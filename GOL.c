#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>

#define ROOT 0
#define ALIVE 1
#define DEAD 0
#define BIGPRIME 68111 //some reasonably big prime number

// Joshua Feltman
// Matthew Johnson
// CS411 Project 2

//commonly referenced integers, changed from main to globals for simplicty
int width, genCount, displayCount, rank, p, miniMatrixSize;

long initialGenerationTime;
long totalCommunicationTime; //all communication steps
long totalGenerationTime; //G sums of time taken to generate

int main(int argc,char *argv[])
{
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&p);

    srand(time(NULL));
    totalGenerationTime = 0L;
    totalCommunicationTime = 0L;
    initialGenerationTime = 0L;

    printf("my rank=%d\n",rank);
    printf("Rank=%d: number of processes =%d\n",rank,p);
    
    char *tempWidth = argv[1];
    char *tempG = argv[2];

    width = atoi(tempWidth);
    genCount = atoi(tempG);
    printf("Width(N) = %d, G = %d\n", width, genCount);

    displayCount = 1; //we will simply display after each generation, e.g. 2 would be after every other gen

    printf("Checking that the number of processors is power of 2 and N is divisible by p...\n");
    if(p != 0 && (p & (p - 1)) == 0 && width % p == 0) 
    {
        printf("Pass.\n");
    }
    else
    {
        printf("Fail. Aborting.\n");
        return 1;
    }
    
    //nxn matrix. Block composition states that a processor owns [width/p][width] cells (or [width][width/p] cells) of the matrix
    //no need to broadcast p since its global now, begin execution
    int miniMatrix[width/p][width];    
    miniMatrixSize = sizeof(miniMatrix);
    
    GenerateInitialGoL(miniMatrix);
    Simulate(miniMatrix);

    MPI_Barrier(MPI_COMM_WORLD);

    printf("RANK: %d | CLOCK Total Communication Time: %ld usec\n", rank, totalCommunicationTime); 
    printf("RANK: %d | CLOCK Total Computation Time: %ld usec\n", rank, (totalGenerationTime+initialGenerationTime)-totalCommunicationTime);
    printf("RANK: %d | CLOCK Total Run Time: %ld usec\n", rank, totalGenerationTime+initialGenerationTime); 
    printf("RANK: %d | CLOCK Average Generation Time: %ld usec\n", rank, totalGenerationTime/genCount);  
    
    MPI_Finalize();

    return 0;
}

void GenerateInitialGoL(int miniMatrix[][width]){
    //Two Step function
    //Step One: rank 0 generates p random prime numbers and <i>distributes</i> them to the other ranks.
    //Step Two: Each rank locally generates (using the passed seed) a distinct sequence of values to fill in the matrix as alive or dead. 

    //Step One
    int seed = 0; //each processors seed after the random numbers are distributed
    int randomNumbers[p]; //List of random numbers to be sent to each process. Buffer is size p for p processes 
    struct timeval start, end, start1, end1;
    
    gettimeofday(&start, NULL);
    if(rank == ROOT){
        int i = 0; //must define outside for loop in MPI
        for(i; i < p; i++){
            randomNumbers[i] = rand() % BIGPRIME + 1; //rand starts at 0, off by 1
        }
    }
    
    gettimeofday(&start1, NULL);
    // Broadcast i random numbers to the i processors
    MPI_Bcast(randomNumbers, p, MPI_INT, ROOT, MPI_COMM_WORLD);
    gettimeofday(&end1, NULL);
    
    long commStep = ((end1.tv_sec * 1000000 + end1.tv_usec) - (start1.tv_sec * 1000000 + start1.tv_usec));
    totalCommunicationTime += commStep;

    // ith random number handed over to rank i
    seed = randomNumbers[rank];
    printf("RANK: %d recieves random number: %d\n", rank, seed);

    //Step Two
    srand(time(NULL) + seed); //redefine random with the seed since every rank after root was getting the same value 
    
    int *tempMini = miniMatrix[0];
    int index=0;
    // Set tempMini with all the Alive/Dead values then we can just copy it into miniMatrix
    // (width^2)/p is number of elements in miniMatrixi
    for(index; index < (width*width)/p; index++){
        int parity = rand() % BIGPRIME + 1;

        //50-50 chance of a cell being alive or dead based on if the status of even or odd (parity)
        if(parity % 2 == 0){
            tempMini[index] = ALIVE;  
        }
        else{
            tempMini[index] = DEAD;
        }
    }

    // Copy tempMini into miniMatrix
    memcpy(miniMatrix, tempMini, miniMatrixSize);
    
    gettimeofday(&end, NULL);
    initialGenerationTime += ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
    
    //Verification Step, -1 represents the initialization
    DisplayGoL(miniMatrix, -1);
}

// Return the new state of a given cell
// Need to get the cell's 8 neighbors (will need to send/recv rows from other processors since
//  we're using a TORUS topology)
int DetermineState(int miniMatrix[][width], int row, int col) {
    struct timeval start, end;
    
    MPI_Status status;
    
    int tempFirstRow[width], tempLastRow[width];
    int numTempMatrixRows = (width / p) + 2;    
    int tempMatrix[numTempMatrixRows][width];
    // Need to send 2 miniMatrix rows (first and last) to both p-1 and p+1
    // Need to recieve 2 rows from both p-1 and p+1
    
    // Since we are using a TORUS topology, the ranks need wrap around e.g. p=4 then rank0 - 1 = rank3 and rank3 + 1 = rank0
    // REF: https://codereview.stackexchange.com/questions/57923/index-into-array-as-if-it-is-circular
    int prevRank = (rank - 1 % p + p) % p;
    int nextRank = (rank + 1) % p;

    gettimeofday(&start, NULL);
    MPI_Send(miniMatrix[0], width, MPI_INT, prevRank, 0, MPI_COMM_WORLD); // send first row
    MPI_Recv(tempLastRow, width, MPI_INT, nextRank, 0, MPI_COMM_WORLD, &status);
   
    MPI_Send(miniMatrix[(width/p) - 1], width, MPI_INT, nextRank, 0, MPI_COMM_WORLD); // send last row
    MPI_Recv(tempFirstRow, width, MPI_INT, prevRank, 0, MPI_COMM_WORLD, &status);
    gettimeofday(&end, NULL);
   
    long commStep = ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
    totalCommunicationTime += commStep;
    
    // Construct tempMatrix, copy in new first and last row
    memcpy(tempMatrix[0], tempFirstRow, sizeof(tempFirstRow));
    memcpy(tempMatrix[1], miniMatrix, miniMatrixSize);
    memcpy(tempMatrix[numTempMatrixRows - 1], tempLastRow, sizeof(tempLastRow));

    // Debugging helper to make sure tempMatrix was assembled correctly
//    MPI_Barrier(MPI_COMM_WORLD);
//    int x, y;
//    printf("RANK: %d\n", rank);
//    for(x = 0; x < (width/p) + 2; x++) {
//        for(y = 0; y < width; y++) {
//            printf("%d ", tempMatrix[x][y]);
//        }
//        printf("\n");
//    }
  
    row += 1; // new row in tempMatrix, shift by 1 since we added a row at the beginning
    int aliveCount = 0;
    
    // Same as when we send/recv from ranks
    int prevCol = (col - 1 % width + width) % width;
    int nextCol = (col + 1) % width;

    // Add up a cells 8 neighbors
    aliveCount += tempMatrix[row - 1][prevCol];
    aliveCount += tempMatrix[row - 1][col];
    aliveCount += tempMatrix[row - 1][nextCol];
        
    aliveCount += tempMatrix[row + 1][prevCol];
    aliveCount += tempMatrix[row + 1][col];
    aliveCount += tempMatrix[row + 1][nextCol];
        
    aliveCount += tempMatrix[row][prevCol];
    aliveCount += tempMatrix[row][nextCol];
    
    //printf("RANK: %d [%d][%d] ALIVECOUNT = %d\n", rank, row-1, col, aliveCount);

    if (aliveCount >= 3 && aliveCount <= 5) {
        return ALIVE;
    } 
    else {
        return DEAD;
    }
}

void Simulate(int miniMatrix[][width]) {
    int i, j;
    int tempMiniMatrix[width/p][width];
    struct timeval start, end, communicationStart, communicationEnd;

    for(i = 0; i < genCount; i++) {
        // Start timer for total runtime of simulate, cant do outside for loop since we dont time display
        gettimeofday(&start, NULL);

        gettimeofday(&communicationStart, NULL);        
        // Make sure each process finishes the generation before starting the next one
        MPI_Barrier(MPI_COMM_WORLD);
        gettimeofday(&communicationEnd, NULL);
    
        long commStep = ((communicationEnd.tv_sec * 1000000 + communicationEnd.tv_usec) - (communicationStart.tv_sec * 1000000 + communicationStart.tv_usec));
        totalCommunicationTime += commStep;

        // (width^2)/p is number of elements in the miniMatrix
        for(j = 0; j < (width * width)/p; j++) {
            int curRow = j / width;
            int curCol = j % width;

            int newState = DetermineState(miniMatrix, curRow, curCol);

            // Make a copy of the MiniMatrix with the new states
            tempMiniMatrix[curRow][curCol] = newState;
        }
        
        // After check all the states, update the miniMatrix with the new states
        // This is because we use the prev gen for all of the checking so we cant update the board as we 
        // check or it will corrupt the game
        memcpy(miniMatrix, tempMiniMatrix, sizeof(tempMiniMatrix));

        gettimeofday(&end, NULL);
        totalGenerationTime += ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));
       
        // Display buffer, dont record its time. 
        if(i % displayCount == 0) {
            DisplayGoL(miniMatrix, i);
        } 
    }
}

void DisplayGoL(int miniMatrix[][width], int generation){
    int fullMatrix[width][width];

    MPI_Barrier(MPI_COMM_WORLD);

    // Printing out memory layout of miniMatrix for debugging help
//    int *x = miniMatrix[0];
//    int i = 0;
//    printf("RANK: %d - ", rank);
//    for(i; i < (width*width/p); i++) {
//        printf("%d ", x[i]);
//    }
//    printf("\n");

    MPI_Gather(miniMatrix[0], (width * width)/p, MPI_INT, fullMatrix[0], (width * width)/p, MPI_INT, ROOT, MPI_COMM_WORLD);

    // Only root displays the matrix since it gathered every other miniMatrix
    if (rank == ROOT) {
        outputMatrix(fullMatrix, generation);
    }
}

void outputMatrix(int fullMatrix[][width], int generation) {
    int i, j;
    
    printf("===== Matrix at Generation %d =====\n", generation + 1);
    for(i = 0; i < width; i++) {
        for(j = 0; j < width; j++) {
            if (fullMatrix[i][j] == ALIVE) {
                printf("A ");
            } else {
                printf("D ");
            }
        }
        printf("\n");
    }
    printf("==================================\n");
}

