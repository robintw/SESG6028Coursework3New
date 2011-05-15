
/* 

Code that solves Laplace's equation in 3D for Dirichlet boundary conditions set in
grid_set_boundary.

Note this only work for NON-periodic boundary conditions

Written by I.J.Bush March 2011 - if it breaks your computer or aught else, tough.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "grid.h"
#include "array_alloc.h"
#include "timer.h"
#include <mpi.h>

int grid_init( int ng[ 3 ], struct grid *g )
{

  /* Initialise the grid 
     Arguments are:
     ng        - The total number of points in the grid
     grid      - The derived type holding all the data about the grid

  */

  int i;
  int npx = 3;
  int npy = 2;
  int npz = 1;

  int rank;
  int periods[3];
  int dim_size[3];
  int coords[3];

  MPI_Comm cart_comm;

  /* Store the size of the whole grid */
  for( i = 0; i < 3; i++ ){
    g->whole_size[ i ] = ng[ i ];
  }
  
  dim_size[0] = npz;
  dim_size[1] = npy;
  dim_size[2] = npx;
  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Create the communicator */	
  MPI_Cart_create(MPI_COMM_WORLD, 3, dim_size, periods, 1, &cart_comm);
	
  /* Get our co-ordinates within that communicator */
  MPI_Cart_coords(cart_comm, rank, 3, coords);
    
  /* Calculate how large a chunk we've got here - this is the size of the chunk that we actually
  want to be able use - so nux, nuy and nuz*/
  g->nuz = ceil(g->whole_size[0] / (float) npz);
  g->nuy = ceil(g->whole_size[1] / (float) npy);
  g->nux = ceil(g->whole_size[2] / (float) npx);
  
  if (coords[0] == (npz - 1))
  {
		/* We're at the far end of z */
		g->nuz = g->whole_size[0] - (g->nuz * (npz - 1));
  }
  if (coords[1] == (npy - 1))
  {
		/* We're at the far end of y */
		g->nuy = g->whole_size[1] - (g->nuy * (npy - 1));
  }
  if (coords[2] == (npx - 1))
  {
		/* We're at the far end of x */
		g->nux = g->whole_size[2] - (g->nux * (npx - 1));
  }
  
  g->pz = coords[0];
  g->py = coords[1];
  g->px = coords[2];
  
  g->npz = npz;
  g->npy = npy;
  g->npx = npx;
  
  /* We want each array to actually have two extra elements in each direction so we can
  store halos or boundary conditions at each end. */
  g->nx = g->nux + 2;
  g->ny = g->nuy + 2;
  g->nz = g->nuz + 2;

  /* Allocate the grid. Note we will need two versions of the data on the grid, one to hold
     the current values, and one to write the results into when we are updating the grid. We swap between
     the two versions as we go from iteration to iteration, the current member of the derived type
     indicating which version is the most up to date. */
  for( i = 0; i < 2; i++ ){

    g->data[ i ] = alloc_3d_double( g->nz, g->ny, g->nx ); 
    if( g->data[ i ] == NULL )
      return EXIT_FAILURE;
  }
  
  printf("Allocated an array of (%d, %d, %d) with usable dimensions of (%d, %d, %d)\n", g->nz, g->ny, g->nx, g->nuz, g->nuy, g->nux);

  /* Which version of the grid is the "current" version. The other we will write
     the next result into */
  g->current = 0;

  /* Zero the timer and the iteration counter*/
  g->t_iter = 0.0;
  g->n_iter = 0;

  return EXIT_SUCCESS;

}

void grid_finalize( struct grid *g ) 
{

  /* Free the grid */

  free_3d_double( g->data[ 1 ], g->nz);
  free_3d_double( g->data[ 0 ], g->nz);

}

void grid_initial_guess( struct grid *g )
{

  /* Initial guess at the solution. We'll use a very simple one ... */

  int i, j, k, l;

    for( i = 0; i < 2; i++ ){
      for( j = 0; j < g->nz; j++ ){
	for( k = 0; k < g->ny; k++ ){
	  for( l = 0; l < g->nx; l++ ){
	    g->data[ i ][ j ][ k ][ l ] = 0.0;
	  }
	}
      } 
    }

}

void grid_set_boundary( struct grid *g )
{
  /* Set the boundary conditions for the grid. We're using Dirichlet
     boundary conditions, so we need to set every point on the edge of the grid.

     As a simple example set all the faces to zero except the bottom xy face, which we set 
     to unity */

  int x, y, z, version;

  /* Set each face of the cuboid in turn */
  /* Also remember that we need to do do it for both versions */
  
  if (g->px == 0)
  {
  	/* We're at the LOW X face of the grid. So set it to zero. */
  	for (version = 0; version < 2; version++)
  	{
  		for (z = 0; z < g->nz; z++)
  		{
  			for (y = 0; y < g->ny; y++)
  			{
  				g->data[version][z][y][0] = 0.0;
  			}
  		}
  	}
  }
  
  if (g->px == (g->npx - 1))
  {
  	/* We're at the HIGH X face of the grid. So set it to zero */
  	for (version = 0; version < 2; version++)
  	{
  		for (z = 0; z < g->nz; z++)
  		{
  			for (y = 0; y < g->ny; y++)
  			{
  				g->data[version][z][y][g->nx-1] = 0.0;
  			}
  		}
  	}
  }
  
  if (g->py == 0)
  {
  	/* We're at the LOW Y face of the grid. So set it to zero */
  	for (version = 0; version < 2; version++)
  	{
  		for (z = 0; z < g->nz; z++)
  		{
  			for (x = 0; x < g->nx; x++)
  			{
  				g->data[version][z][0][x] = 0.0;
  			}
  		}
  	}
  }
  
  if (g->py == (g->npy - 1))
  {
  	/* We're at the HIGH Y face of the grid. So set it to zero */
  	for (version = 0; version < 2; version++)
  	{
  		for (z = 0; z < g->nz; z++)
  		{
  			for (x = 0; x < g->nx; x++)
  			{
  				g->data[version][z][g->ny-1][x] = 0.0;
  			}
  		}
  	}
  }
  
  if (g->pz == 0)
  {
  	for (version = 0; version < 2; version++)
  	{
  		for (y = 0; y < g->ny; y++)
  		{
  			for (x = 0; x < g->nx; x++)
  			{
  				g->data[version][0][y][x] = 0.0;
  			}
  		}
  	}
  }
  
  if (g->pz == (g->npz - 1))
  {
  	for (version = 0; version < 2; version++)
  	{
  		for (y = 0; y < g->ny; y++)
  		{
  			for (x = 0; x < g->nx; x++)
  			{
  				g->data[version][g->nz-1][y][x] = 0.0;
  			}
  		}
  	}
  }
}

double grid_update( struct grid *g ){

  /* Perform the grid update */

  #define ONE_SIXTH 0.166666666666666666666666666666666666666666666

  double dg, diff;
  double start, finish;

  int current, update;
  int lb0, lb1, lb2, ub0, ub1, ub2;
  int i, j, k;

  /* Work out which version of the grid hold the current values, and
     which we will write the update into */
  switch( g->current ) {
  case 0:  
    current = 0;
    update  = 1; 
    break;  
  case 1:  
    current = 1;
    update  = 0; 
    break;
  default: 
    fprintf( stderr, "Internal error: impossible value for g->current\n" );
    exit( EXIT_FAILURE );
  }

  /* Bounds for the loops. Remember we should not update the
     boundaries as they are set by the boundary conditions */
  lb0 = 1;
  lb1 = 1;
  lb2 = 1;

  ub0 = g->nz - 2;
  ub1 = g->ny - 2;
  ub2 = g->nx - 2;

  /* Perform the update and check for convergence  */
  start = timer();
  dg = 0.0;
  for( i = lb0; i <= ub0; i++ ) {
    for( j = lb1; j <= ub1; j++ ) {
      for( k = lb2; k <= ub2; k++ ) {
	g->data[ update ][ i ][ j ][ k ] = 
	  ONE_SIXTH * ( g->data[ current ][ i + 1 ][ j     ][ k     ] +
			g->data[ current ][ i - 1 ][ j     ][ k     ] +
			g->data[ current ][ i     ][ j + 1 ][ k     ] +
			g->data[ current ][ i     ][ j - 1 ][ k     ] +
			g->data[ current ][ i     ][ j     ][ k + 1 ] +
			g->data[ current ][ i     ][ j     ][ k - 1 ] );
	diff = fabs( g->data[ update ][ i ][ j ][ k ] - g->data[ current ][ i ][ j ][ k ] );
	dg = dg > diff ? dg : diff;
      }
    }
  }
  finish = timer();
  g->t_iter += finish - start;

  /* Update the iteration counter */
  g->n_iter++;

  /* The updated grid is now the current grid, so swap over */
  g->current = update;

  return dg;

}

double grid_checksum( struct grid g ){

  /* Simple checksum over the grid to help checking solution
     i..e simply add up all the grid points */

  double sum;

  int i, j, k;

  sum = 0.0;
  for( i = 0; i < g.nz; i++ ) {
    for( j = 0; j < g.ny; j++ ) {
      for( k = 0; k < g.nx; k++ ) {
	sum += g.data[ g.current ][ i ][ j ][ k ];
      }
    }
  }

  return sum;

}

void grid_print_times( struct grid g ){

  /* Report the measured timing data */

  fprintf( stdout, "Timing breakdown for the grid operations:\n" );
  fprintf( stdout, "                                  Total     Per iteration\n" );
  fprintf( stdout, "Iteration time        :      %12.4f   %12.8f\n", g.t_iter, g.t_iter / g.n_iter );
  fprintf( stdout, "\n" );

}
