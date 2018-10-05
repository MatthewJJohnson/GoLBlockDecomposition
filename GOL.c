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

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    
    printf("Checking that the number of processors is power of 2...\n");
    if(p!=0 && ((p & (p -1)) == 0)
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
    
    int width;
    int genCount;

    if(rank == ROOT)
    {
        userInput(&width, &genCount, p);
        for (int i=1; i<p; i++)
        {
            MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        printf("Process %d received the width of the board as %d from process 0\n", rank, width);
    }
    GenerateInitialGoL();
    
}

void userInput(int *width, int *genCount)
{
    int pass = 0;
    while(pass==0)
    {
        printf("Enter the desired width/column count (n): ");
        scanf("%d", *width);
        printf("\nChecking that n>p and n is divisible by p...\n");
        if(*width%p == 0)
        {
            pass = 1;
            printf("Pass.\n");
        }
        else
        {
            printf("Fail.\n");
        }
    }
    printf("Enter the number of generations the GOL will have: ");
    scanf("%d", *genCount);
    printf("\n");
}

void GenerateInitialGoL(int *rank, int *){
    
}

void Simulate(){

}

void DisplayGoL(){

}




