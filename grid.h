
/* Include file for grid.c - for more details see there */

/* Written by I.J.Bush March 2011 - if it breaks your computer or aught else, tough. */ 

/* The grid and associated data. Note we will need two versions of the data on the grid, one to hold
   the current values, and one to write the results into when we are updating the grid. Hence the
   grid has 4 dimensions not 3, the last dimension being the two versions of the grid. We swap between
   the two versions as we go from iteration to iteration, the current member of the derived type
   indicating which version is the most up to date.*/
struct grid {
  
  double ***data[ 2 ]; /* The (two versions of) the data on the grid */
  int current;         /* Which version of the grid we are currently using */
  int whole_size[ 3 ];         /* The size of the overall grid */
  
  int nx, ny, nz; /* The size of the chunk given to this process */
  int nux, nuy, nuz; /* The actual usable dimensions of this chunk - ignoring the edge bits used for halos/boundaries */
  
  int px, py, pz; /* Positions in the grid */
  int npx, npy, npz; /* Number of processors in the grid */
  
  int north, south, east, west, up, down; /* Neighbours */
  
  int n_iter;          /* The number of iterations */
  double t_iter;       /* The time in the iterations */

};

/* The prototypes */
int grid_init( int ng[ 3 ], struct grid *g );    /* Initialise the grid */
void grid_finalize( struct grid *g );            /* Finalise the grid   */ 
void grid_initial_guess( struct grid *g );       /* The initial guess at the solution */
void grid_set_boundary( struct grid *g );        /* Set the bondary conditions */
double grid_update( struct grid *g );            /* Perform one iteration of the grid update */
double grid_checksum( struct grid g );           /* provide a checksum as a check of the answer */
void grid_print_times( struct grid g );          /* Report timing data */

