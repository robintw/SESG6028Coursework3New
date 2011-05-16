
/* Diver routine to solve Laplace's equaiton in 3D */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "grid.h"
#include "timer.h"
#include <mpi.h>

int main(int argc, char **argv){

  struct grid g;

  time_t tp;

  double tol, dg, max_dg;
  double check, total_check;
  double start = 0.0;
  double finish = 0.0;

  int ng[ 3 ];

  int max_iter;
  int error;
  int converged;
  int retval;
  int i;
  
  int npx, npy, npz;
  
  int rank, nprocs;

  /* Ever the optimist ... Assume the code has failed until we see otherwise */
  retval = EXIT_FAILURE;
  
  /* Initialise the MPI subsystem */
  MPI_Init(&argc, &argv);
  
  /* Get the number of processes we have access to and the rank of the current process */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
	/* Read in the input values */
	fprintf( stdout, "How big should the grid be (nx,ny,nz)?\n" );
	if( fscanf ( stdin , "%i %i %i", &ng[ 0 ], &ng[ 1 ], &ng[ 2 ] ) != 3 ) {
	fprintf( stderr, "ERROR: Failed to read input: Grid points line incorrect\n" );
	exit( retval );
	}
	fprintf( stdout, "What tolerance should I converge to?\n" );
	if( fscanf ( stdin , "%lf", &tol ) != 1 ) {
	fprintf( stderr, "ERROR: Failed to read input: Tolerance line incorrect\n" );
	exit( retval );
	}
	fprintf( stdout, "What is the maximum number of iterations?\n" );
	if( fscanf ( stdin , "%i", &max_iter ) != 1 ) {
	fprintf( stderr, "ERROR: Failed to read input: Number of iterations line incorrect\n" );
	exit( retval );
	}
	
	fprintf( stdout, "What is the number of processes in the X direction\n" );
	if( fscanf ( stdin , "%i", &npx ) != 1 ) {
		fprintf( stderr, "ERROR: Failed to read input: Number of processes in X direction line incorrect\n" );
		exit( retval );
	}
	
	fprintf( stdout, "What is the number of processes in the Y direction\n" );
	if( fscanf ( stdin , "%i", &npy ) != 1 ) {
		fprintf( stderr, "ERROR: Failed to read input: Number of processes in Y direction line incorrect\n" );
		exit( retval );
	}
	
	fprintf( stdout, "What is the number of processes in the Z direction\n" );
	if( fscanf ( stdin , "%i", &npz ) != 1 ) {
		fprintf( stderr, "ERROR: Failed to read input: Number of processes in Z direction line incorrect\n" );
		exit( retval );
	}
	
	fprintf( stdout, "\n" );
	
	/* Report the input data */
	fprintf( stdout, "The calculation is for a %i * %i * %i grid\n"
	   "The convergance tolerance is %f\n"
	   "The maximum number of iterations is %i\n", 
		   ng[ 0 ], ng[ 1 ], ng[ 2 ], tol, max_iter );
	tp = time( &tp );
		fprintf( stdout, "Running on %d processors\n", nprocs);
	fprintf( stdout, "Running this job at %s\n", ctime( &tp ) );
  }

  /* Broadcast the values we've just got in from the user to the other processes */
  MPI_Bcast(&ng, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max_iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&npx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npy, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npz, 1, MPI_INT, 0, MPI_COMM_WORLD);



  /* Initialise the grid */
  error = grid_init( ng, npx, npy, npz, &g );

  /* Set the inital guess at the solution */
  grid_initial_guess( &g );

  /* Set the boundary conditions */
  grid_set_boundary( &g );

  /* Loop updating the grid until the change is sufficently small to consider
     the calculation converged */
       
  /* Start the timer */
  if (rank == 0)
  {
   	start = timer();
  }

  /* Assume we haven't converged */
  converged = 0;
  
  
  for( i = 1; i <= max_iter; i++ )
  {
  	/* Do the actual update */
	dg = grid_update( &g );

   	/* Pass the maximum change (returned from grid_update) back to the root process
	so that it can then check to see if it is smaller than the tolerance */
	MPI_Reduce(&dg, &max_dg, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	/* On the root process ONLY, decide whether we've converged */
	if (rank == 0)
	{
		/* Let the user know how it's going */
		fprintf( stdout, "Iter %5i Max change %20.12f\n", i, max_dg );
		
		/* If the reduced difference is less than the tolerance then say that we've converged */
		if (max_dg < tol)
		{
			converged = 1;
		}
	}
	
	/* Send a value to all of the processes to say if we've converged or not */
	MPI_Bcast(&converged, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	/* If we have, then exit the loop (all processes will do this) */
	if (converged == 1)
	{
		break;
	}
   }
    
	/* Get the time now we've finished updating the grid */
   if (rank == 0)
   {
   		finish = timer();
   }

    /* Add up all the grid points - can be used as a simple
     check that things have worked */
    
    /* Get the checksum */
    check = grid_checksum( g );
    
    total_check = 0;
    /* Sum all of the sub-checksums by using a MPI_Reduce */
    MPI_Reduce(&check, &total_check, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	/* Print out the results */
	if (rank == 0)
	{
		/* Report what happened */
		fprintf( stdout, "\n" );
		if( converged ) {
		  fprintf( stdout, "The calculation converged\n" );
		}
		else{
		  fprintf( stdout, "The calculation FAILED to converged\n" );
		}
		fprintf( stdout, "\n" );
		fprintf( stdout, "The calculation took %f seconds\n", finish - start );
		fprintf( stdout, "The check sum of the grid is %20.12g\n", total_check );
		fprintf( stdout, "\n" );
		grid_print_times( g );
	
		/* Calculation worked so tell the world */
		retval = EXIT_SUCCESS;
	}

  /* finalise the grid */
  grid_finalize( &g );

  MPI_Finalize();
  return retval;

}
