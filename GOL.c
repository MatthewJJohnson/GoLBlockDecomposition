#include <stdio.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>
#include <math.h>

#define ROOT 0

// Joshua Feltman
// Matthew Johnson
// CS411 Project 2

int main(int argc,char *argv[])
{
    int rank,p;
    int width, genCount;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    
    printf("Checking that the number of processors is power of 2...\n");
    if(p != 0 && (p & (p - 1)) == 0) 
    {
        printf("Pass.\n");
    }
    else
    {
        printf("Fail.\n");
        return 1;
    }

    printf("my rank=%d\n",rank);
    printf("Rank=%d: number of processes =%d\n",rank,p);
    
    char *tempWidth = argv[1];
    char *tempG = argv[2];

    width = atoi(tempWidth);
    genCount = atoi(tempG);

    printf("N = %d, G = %d\n", width, genCount);

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
    
}

void Simulate(){

}

void DisplayGoL(){

}

