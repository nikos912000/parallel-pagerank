/********************************************************************/
/*    Pagerank project 2014 - Serial version                        */
/*    	*based on Cleve Moler's matlab implementation               */
/*                                                                  */
/*    Implemented by Nikos Katirtzis (nikos912000)                  */
/********************************************************************/

/******************** Includes - Defines ****************/
#include "pagerank_serial.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>

/******************** Defines ****************/
// Number of nodes
int N;

// Convergence threashold and algorithm's parameter d  
double threshold, d;

// Table of node's data
Node *Nodes;


/***** Read graph connections from txt file *****/	

void Read_from_txt_file(char* filename)
{
    
    FILE *fid;

    int from_idx, to_idx;
	int temp_size;

    fid = fopen(filename, "r");
   	if (fid == NULL){printf("Error opening data file\n");}

	while (!feof(fid))
	{
		
		if (fscanf(fid,"%d\t%d\n", &from_idx,&to_idx))
		{
			Nodes[from_idx].con_size++;
			temp_size = Nodes[from_idx].con_size;
			Nodes[from_idx].To_id =(int*) realloc(Nodes[from_idx].To_id, temp_size * sizeof(int));
			Nodes[from_idx].To_id[temp_size - 1] = to_idx;
		}
	}

	printf("End of connections insertion!\n");

	fclose(fid);

}

/***** Read P vector from txt file*****/	

void Read_P_from_txt_file()
{

	FILE *fid;
	double temp_P;
	int index = 0;

    fid = fopen("P.txt", "r");
   	if (fid == NULL){printf("Error opening the Probabilities file\n");}

	while (!feof(fid))
	{
		// P's values are double!
		if (fscanf(fid,"%lf\n", &temp_P))
		{
			Nodes[index].p_t1 = temp_P;
			index++;	   
		}
	}
	printf("End of P insertion!");

	fclose(fid);	

}


/***** Read E vector from txt file*****/	

void Read_E_from_txt_file()
{

	FILE *fid;
	double temp_E;
	int index = 0;

    fid = fopen("E.txt", "r");
   	if (fid == NULL){printf("Error opening the E file\n");}

	while (!feof(fid))
	{
		// E's values are double!
		if (fscanf(fid,"%lf\n", &temp_E))
		{
			Nodes[index].e = temp_E;
			index++;   
		}
	}
	printf("End of E insertion!");

	fclose(fid);	

}

/***** Create P and E with equal probability *****/

void Random_P_E()
{

   	int i;
    // Sum of P (it must be =1)
    double sum_P_1 = 0;
    // Sum of E (it must be =1)
    double sum_E_1 = 0; 
       
    // Arrays initialization
    for (i = 0; i < N; i++)
    {
        Nodes[i].p_t0 = 0;
        Nodes[i].p_t1 = 1;
        Nodes[i].p_t1 = (double) Nodes[i].p_t1 / N;

        sum_P_1 = sum_P_1 + Nodes[i].p_t1;
        
		Nodes[i].e = 1;
        Nodes[i].e = (double) Nodes[i].e / N;
        sum_E_1 = sum_E_1 + Nodes[i].e;
    }

    // Assert sum of probabilities is =1
    
    // Print sum of P (it must be =1)
    //printf("Sum of P = %f\n",sum_P_1);
    
    // Exit if sum of P is !=1
    assert(sum_P_1 = 1);
    
    //printf("\n");
    
    // Print sum of E (it must be =1)
    //printf("Sum of E = %f\n",sum_E_1);
    
    // Exit if sum of Pt0 is !=1
    assert(sum_E_1 = 1);

}


/***** Main function *****/   

int main(int argc, char** argv)
{
	// Check input arguments
	if (argc < 5)
	{
		printf("Error in arguments! Three arguments required: graph filename, N, threshold and d\n");
		return 0;
	}

	// get arguments 
	char filename[256];
	strcpy(filename, argv[1]);
	N = atoi(argv[2]);
	threshold = atof(argv[3]);
	d = atof(argv[4]);
    
    int i, j, k;
	double totaltime;

	// a constant value contributed of all nodes with connectivity = 0
	// it's going to be addes to all node's new probability
	double sum = 0;
    
	// Allocate memory for N nodes
    Nodes = (Node*) malloc(N * sizeof(Node));
    
    for (i = 0; i < N; i++)
	{
		Nodes[i].con_size = 0;
        Nodes[i].To_id = (int*) malloc(sizeof(int));
    }

    Read_from_txt_file(filename);
    
	// set random probabilities
    Random_P_E();
	// OR read probabilities from files
	//Read_P_from_txt_file();
	//Read_E_from_txt_file();
    
	gettimeofday(&start, NULL);
    
	
    /********** Start of algorithm **********/
    
    // Iterations counter
    int iterations = 0;
    int index;
    
    // Or any value > threshold
    double max_error = 1;
    
    printf("\nSerial version of Pagerank\n");

    // Continue if we don't have convergence yet
    while (max_error > threshold)
    {
    	sum = 0;
		
		// Initialize P(t) and P(t + 1) values
        for (i = 0; i < N; i++)
        {
            // Update the "old" P table with the new one 
            Nodes[i].p_t0 = Nodes[i].p_t1;  
            Nodes[i].p_t1 = 0;
        }

        // Find P for each webpage
        for (i = 0; i < N; i++)
        {
            
            if (Nodes[i].con_size != 0)
            {
                
                // Compute the total probability, contributed by node's neighbors
                for (j = 0; j < Nodes[i].con_size; j++)
                {
					index = Nodes[i].To_id[j];	
					Nodes[index].p_t1 = Nodes[index].p_t1 + (double) Nodes[i].p_t0 / Nodes[i].con_size;
				}

            }
            
            else
            {
				// Contribute to all
                sum = sum + (double)Nodes[i].p_t0 / N;
            }
            
        }
        
		max_error = -1;
        
		// Compute the new probabilities and find maximum error
        for (i = 0;i < N; i++)
        {
            Nodes[i].p_t1 = d * (Nodes[i].p_t1 + sum) + (1 - d) * Nodes[i].e;
            
            if (fabs(Nodes[i].p_t1 - Nodes[i].p_t0) > max_error)
            {
                max_error = fabs(Nodes[i].p_t1 - Nodes[i].p_t0);
            }
            
        }
        
        printf("Max Error in iteration %d = %f\n", iterations+1, max_error);
        iterations++;
    }

	gettimeofday(&end, NULL);
    
	printf("\n");
  
	// Print final probabilitities
    /*for (i = 0; i < N; i++)
    {
        printf("P_t1[%d] = %f\n",i,Nodes[i].p_t1);
    }
    
    printf("\n");*/
    
    // Print no of iterations
    printf("Total iterations: %d\n", iterations);
	
    totaltime = (((end.tv_usec - start.tv_usec) / 1.0e6 + end.tv_sec - start.tv_sec) * 1000) / 1000;

	printf("\nTotaltime = %f seconds\n", totaltime);
    
    printf("End of program!\n");
    
    return (EXIT_SUCCESS);
}