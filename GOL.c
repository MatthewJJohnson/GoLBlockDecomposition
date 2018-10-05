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
    
    //no need to broadcast p since its global now, begin execution
    int miniMatrix[width][width/p];    
    GenerateInitialGoL(miniMatrix);

    MPI_Finalize();    
}

void GenerateInitialGoL(int miniMatrix[n][]){
//Two Step function
//Step One: rank 0 generates p random prime numbers and <i>distributes</i> them to the other ranks. Instuctor hinted at the MPI function, 'MPI_Scatter'
//Step Two: Each rank locally generates (using the passed seed) a distinct sequence of values to fill in the matrix as alive or dead. 

//Step One
    //REF: www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node72.html
    int scatteredSeed;//save our seed from the scatter
    int toBeScattered[p];//MPI_Scatter takes a sendbuf (toBeScattered) and splits the buffer between the processes. buffer is size p for p processes 
    
    srand(time(NULL));
    if(rank == ROOT){
        int i=0;//must define outside for loop in MPI
        for(i; i<p;i++){
            toBeScattered[i] = rand() % BIGPRIME+1; //rand starts at 0, off by 1

        }
    }
    MPI_Scatter(toBeScattered, 1, MPI_INT, &scatteredSeed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
//Step Two
    srand(time(NULL)+scatteredSeed);//redefine random seed
    int index=0;
    int parity;
    for(index;index<(width/p)*width;index++){
        parity = rand()%BIGPRIME+1;
        //50-50 chance of a cell being alive or dead based on if the status of even or odd (parity)
        printf("Processor %d generated Cell[%d][%d] as ", rank, index%width, index/width);
        if(parity % 2 != 0){
            miniMatrix[index%width][index/width] = ALIVE;
            printf("ALIVE\n");
        }
        else{
            miniMatrix[index%width][index/width] = DEAD;
            printf("DEAD\n");
        }
    }

//Verification Step
    //DisplayGoL();
}

void Simulate(){

}

void DisplayGoL(){

}

