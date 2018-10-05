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
    //we will implement the [width/p][width] method. (row decomposition)

    int miniMatrix[width/p][width];

    if(rank == ROOT)
    {
	    int i = 1;
        for (i; i < p; i++)
        {
            MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Status status;

        MPI_Recv(&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        printf("Process %d received the width of the board as %d from process 0\n", rank, width);
    }
//    GenerateInitialGoL();

    MPI_Finalize();    
}

void GenerateInitialGoL(int *rank, int *p){
//Two Step function
//Step One: rank 0 generates p random prime numbers and <i>distributes</i> them to the other ranks. Instuctor hinted at the MPI function, 'MPI_Scatter'
//Step Two: Each rank locally generates (using the passed seed) a distinct sequence of values to fill in the matrix as alive or dead. 

//Step One
    
}

void Simulate(){

}

void DisplayGoL(){

}

