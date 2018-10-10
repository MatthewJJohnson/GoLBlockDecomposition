#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>
#include <math.h>

#define ROOT 0
#define ALIVE 1
#define DEAD 0
#define BIGPRIME 68111 //some reasonably big prime number

// Joshua Feltman
// Matthew Johnson
// CS411 Project 2

//commonly referenced integers, changed from main to globals for simplicty
int width, genCount, displayCount, liveCount, deadCount, rank, p;

int main(int argc,char *argv[])
{
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&p);

    srand(time(NULL));

    printf("my rank=%d\n",rank);
    printf("Rank=%d: number of processes =%d\n",rank,p);
    
    char *tempWidth = argv[1];
    char *tempG = argv[2];

    width = atoi(tempWidth);
    genCount = atoi(tempG);
    displayCount = genCount;//we will simply display after each generation
    printf("N = %d, G = %d\n", width, genCount);

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
    //we will implement the [width][width/p] method. (col decomposition)
    // ** See note above DisplayGoL, switching to row decomposition **

    //no need to broadcast p since its global now, begin execution
    int miniMatrix[width/p][width];    
    GenerateInitialGoL(miniMatrix);

    MPI_Finalize();

    return 0;
}

void GenerateInitialGoL(int miniMatrix[][width]){
    //Two Step function
    //Step One: rank 0 generates p random prime numbers and <i>distributes</i> them to the other ranks. Instuctor hinted at the MPI function, 'MPI_Scatter'
    //Step Two: Each rank locally generates (using the passed seed) a distinct sequence of values to fill in the matrix as alive or dead. 

    //Step One
    //REF: www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node72.html
    int scatteredSeed = 0;//save our seed from the scatter
    int toBeScattered[p];//MPI_Scatter takes a sendbuf (toBeScattered) and splits the buffer between the processes. buffer is size p for p processes 
    
    if(rank == ROOT){
        int i = 0; //must define outside for loop in MPI
        for(i; i < p; i++){
            toBeScattered[i] = rand() % BIGPRIME + 1; //rand starts at 0, off by 1
        }
    }
    MPI_Scatter(toBeScattered, 1, MPI_INT, &scatteredSeed, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    printf("RANK: %d recieves random number: %d\n", rank, scatteredSeed);

    //Step Two
    srand(time(NULL) + scatteredSeed); //redefine random with the seed since every rank after root was getting the same value 
    int index=0;
    for(index; index < (width*width)/p; index++){
        int parity = rand() % BIGPRIME + 1;

        //50-50 chance of a cell being alive or dead based on if the status of even or odd (parity)
        printf("Processor %d generated Cell[%d][%d] as ", rank, index/width, index%width);
        if(parity % 2 == 0){
            miniMatrix[index/width][index%width] = ALIVE;
            printf("ALIVE\n");
        }
        else{
            miniMatrix[index/width][index%width] = DEAD;
            printf("DEAD\n");
        }
    }

    //Verification Step
    DisplayGoL(miniMatrix);
}

void Simulate(){

}

// in progress
// *Switching our decomposition from column to row since the memory layout of our mini matrix makes it easier to gather in the 
// full matrix. 
// E.g. miniMatrix[0] = { 1, 2, 3, 4, 5, 6, 7, 8 } can represent a 4x2 Matrix and 2x4 Matrix
//  With a 4x2 Matrix it goes 1 2   and a 2x4 Matrix is 1 2 3 4 
//                            3 4                       5 6 7 8
//                            5 6
//                            7 8
//  When we add these submatrices into a full 4x4 matrix, the 2x4 memory layout doesn't need to be reordered but the 4x2 would 
//  need to be reordered
void DisplayGoL(int miniMatrix[][width]){
    int fullMatrix[width][width];
   
    MPI_Barrier(MPI_COMM_WORLD);

    // Printing out in memory layout of mini matrices for debugging help
    // remove before turning in. 
    int *x = miniMatrix[0];
    int i = 0;
    printf("RANK: %d - ", rank);
    for(i; i < (width*width/p); i++) {
        printf("%d ", x[i]);
    }
    printf("\n");

    // Needed to help get the correct order of sub matrices inside the full Matrix
    if(rank == ROOT) {
         MPI_Gather(miniMatrix[0], (width * width)/p, MPI_INT, fullMatrix[0], (width * width)/p, MPI_INT, ROOT, MPI_COMM_WORLD);
    }
    else {
         MPI_Gather(miniMatrix[0], (width * width)/p, MPI_INT, NULL, 0, MPI_INT, ROOT, MPI_COMM_WORLD);
    }
    
    // Can remove this Barrier when we remove the above debugging code
    MPI_Barrier(MPI_COMM_WORLD);

    // Only root displays the matrix since it gathered every other miniMatrix
    if (rank == ROOT) {
        int i, j;
      
        printf("----COMBINED MATRIX----\n");
        for(i = 0; i < width; i++) {
            for(j = 0; j < width; j++) {    
                printf("%d ", fullMatrix[i][j]);
            }
            printf("\n");
        }
        printf("-----------------------\n");
    }
}

