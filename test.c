
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

  double tol, dg, check;
  double start, finish;

  int ng[ 3 ];

  int max_iter;
  int error;
  int converged;
  int retval;
  int i;
  
  int rank, nprocs;

  /* Ever the optimist ... Assume the code has failed until we see otherwise */
  retval = EXIT_FAILURE;
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  printf("Running with %d processes\n", nprocs);

  if (rank == 0)
  {
	/* Read in the input values */
	fprintf( stdout, "How big should the gird be (nx,ny,nz)?\n" );
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
	fprintf( stdout, "\n" );
	
	/* Report the input data */
	fprintf( stdout, "The calculation is for a %i * %i * %i grid\n"
	   "The convergance tolerance is %f\n"
	   "The maximum number of iterations is %i\n", 
		   ng[ 0 ], ng[ 1 ], ng[ 2 ], tol, max_iter );
	tp = time( &tp );
	fprintf( stdout, "Running this job at %s\n", ctime( &tp ) );
  }


  /* Initialise the grid */
  error = grid_init( ng, &g );

    /* Set the inital guess at the solution */
    grid_initial_guess( &g );

    /* Set the boundary conditions */
    grid_set_boundary( &g );

    /* Loop updating the grid until the change is sufficently small to consider
       the calculation converged */
    start = timer();
    converged = 0;
    for( i = 1; i <= max_iter; i++ ){
      dg = grid_update( &g );
      /* Periodic report of the calculation's status */
      /*if( ( i == 1 || i%10 == 0 ) || dg < tol ) {*/
	fprintf( stdout, "Iter %5i Max change %20.12f\n", i, dg );
      /*}*/
      if( dg < tol ) {
	converged = 1;
	break;
      }
	
    }
    finish = timer();

    /* Add up all the grid points - can be used as a simple
     check that things have worked */
    check = grid_checksum( g );

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
		fprintf( stdout, "The check sum of the grid is %20.12g\n", check );
		fprintf( stdout, "\n" );
		grid_print_times( g );
	
		/* finalise the grid */
		grid_finalize( &g );
	
		/* Calculation worked so tell the world */
		retval = EXIT_SUCCESS;
	}

  MPI_Finalize();
  return retval;

}
