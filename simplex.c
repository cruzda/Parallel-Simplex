#include "simplex.h"

/***************************************************************************//**
/** Initialize Simplex
* Creates & Initializes a new simplex object;  Its address is returned.
*
* The function allocates memory for a new simplex object, sets its boundaries,
* dimensions, and initial points.  It then calculates the function value at 
* those points and returns the address of the new simplex.
*
* @param ndims
*    Specifies the number of dimensions in the simplex.  This is constant 
*    through the life of the simplex.
* @param npts
*    Specifies the number of points in the simplex.  This is constant through
*    the life of the simplex
* @param source
*    Argument tells initialize_simplex() the source of the initial points for
*    the simplex:
*        0 indicates unit vectors 
*        1 indicates random points within the soft boundaries
*        2 indicates points are taken from pts parameter
*    Neither option 0 or 2 are guaranteed to be within the boundaries.
*    Choosing option 1 in combination with bounds==NULL will cause the program
*        to quit with a complaint about arguments.
* @param pts
*    When source = 2 this must point to a block of memory of size 
*    >= npts*ndims*sizeof(double) which contains the initial points of the 
*    simplex.  This memory should be allocated with malloc as the simplex 
*    destructor will attempt to free this memory when it is called.  The memory 
*    block can be thought of as a two dimensional array where each "row" begins 
*    at the index [row*ndims]. Each row contains ndims doubles specifing the 
*    location of one point in the simplex.
* @param bounds
*    Points a block of memory of size >= 3*ndims*sizeof(double) which contains
*    information defining two boundaries for the simplex.  See comments in
*    simplex.h for a full description of how to define these boundaries.
*
*    bounds can also be NULL in which case no boundaries are enforced.
* @param func
*    *func() is the cost function that each point of the simplex will be 
*    evaulated at.  Its arguments are the dimension of the simplex and
*    a pointer to the paramaters at which it should be evaluated.
*
* @return address of the new simplex
*******************************************************************************/
struct simplex *initialize_simplex(
    int const ndims, 
    int const npts, 
    int const source,  /* 0=unit vectors; 1=random; 2=pts argument */
    double * const pts,
    double * bounds,
    double (*const func)(int const, double const * const))
{
  struct simplex *simpx;
  simpx = malloc(sizeof(struct simplex));
  simpx->ndims  = ndims;
  simpx->npts   = npts;
  simpx->borders = bounds;

  /* Each "row" of *points contains a point.  
   * "Columns" contain the point's value for each dimension */
  /* Access element (i,j) of the matrix by p[i*ndim + j] */

  simpx->vals   = (double *) malloc (npts  * sizeof (double));
  simpx->psum   = (double *) malloc (ndims * sizeof (double));

  /* Use Calloc here so that change is initialized to 0s */
  int ptIt, dimIt;
  simpx->change = (double *) calloc (ndims,  sizeof (double)); 
  for(dimIt = 0; dimIt < ndims; dimIt++)
  {
    simpx->change[dimIt] = 0;
  }
  
  assert( simpx->vals && simpx->psum && simpx->change );
  if( source < 2 )
  {
    simpx->points = (double *) malloc (npts  * ndims * sizeof (double));
    assert( simpx->points );
  }

  srand( time(NULL) );
  for (ptIt = 0; ptIt < npts; ptIt++)
    {
      for (dimIt = 0; dimIt < ndims; dimIt++)
	{
	  switch(source)
	  {
	    case 0:
	      /* Initialize simpx->points to unit vectors & origin */
	      simpx->points[ptIt * ndims + dimIt] = 
		     (ptIt == (dimIt + 1) ? 1.0 : 0.0);
	      break;
	    case 1:
	      /* Initialize simpx->points to random points within boundaries */
	      if( bounds == NULL)
	      {
		fprintf(stderr, "borders must be set for random simplex generation\n");
	        exit(0);
	      }
	      simpx->points[ptIt*ndims + dimIt] = ((float)rand() / (float) RAND_MAX) *
		  floor(bounds[2*ndims+dimIt] - bounds[1*ndims+dimIt])
		+ bounds[1*ndims+dimIt];
	      break;
	    case 2:
	      simpx->points = pts;
	      break;
	    default:
	      return NULL;
	  }
	}
      simpx->vals[ptIt] = func(ndims, &simpx->points[ ptIt * ndims ]);
    }
  return simpx;
}


/***************************************************************************//**
* Frees the memory allocated in initialize_simplex() and sets pointer to NULL.
*******************************************************************************/
void destroy_simplex(struct simplex *simpx)
{
  free(simpx->points);
  free(simpx->vals);
  free(simpx->psum);
  free(simpx->borders);
  free(simpx->change);

  free(simpx);

  simpx = NULL;

  return;
}


/***************************************************************************//**
* Find lowest, highest, and 2nd highest vals in simplex 
*
* Function iterates through each point in the simplex to determine the indices
* of the highest, 2nd highest and lowest points.  These indices are updated in
* simpx->ihi, simpx->inhi and simpx->ilo respectively.
*
* @param simpx
*     Address of simplex to find the highest, 2nd higest, and lowest points in.
*******************************************************************************/
void find_hilo(struct simplex * const simpx)
{
  /* Initialize ilo, inhi, and ihi based on first two points */
  if (simpx->vals[0] > simpx->vals[1])
  {
    simpx->ilo  = 1;
    simpx->inhi = 1;
    simpx->ihi  = 0;
  } else {
    simpx->inhi = 0;
    simpx->ilo  = 0;
    simpx->ihi  = 1;
  }

  /* Step through remaining points, if any, updating indexes as needed */
  int ptIt;
  for( ptIt = 0; ptIt < simpx->npts; ptIt++ )
  {
    if ( simpx->vals[ptIt] <= simpx->vals[simpx->ilo] )
    { /* Current point is lower than all previous */
      simpx->ilo = ptIt;
    }
    else if ( simpx->vals[ptIt] > simpx->vals[simpx->ihi] )
    { /* Current point is higher than all previous */
      simpx->inhi = simpx->ihi;  /* Save previous high as next highest */
      simpx->ihi  = ptIt;         /* Update highest point in simplex */
    }
    else if ( simpx->vals[ptIt] > simpx->vals[simpx->inhi] && 
	      ptIt != simpx->ihi )
    { /* Current point is in the set (inhi, ihi) */
      simpx->inhi = ptIt;
    }
  }

  return;
}

/***************************************************************************//**
* Calculate psum for simplex
* Updates psum to contain the sum of each component of each vertex.
*
* @param simpx
*    Contains the address of the simplex whose psum is to be updated.
*******************************************************************************/
void get_psum(struct simplex * const simpx )
{
  int dimIt;
  int ptIt;

  for (dimIt = 0; dimIt < simpx->ndims; dimIt++)
  {
    simpx->psum[dimIt] = 0.0;
    for (ptIt=0; ptIt < simpx->npts; ptIt++)
    {
      simpx->psum[dimIt] += simpx->points[ptIt * simpx->ndims + dimIt];
    }
  }

  return;
}



/***************************************************************************//**
 * Prints the points of the simplex to stdout
*******************************************************************************/
void print_simplex(struct simplex const * const simpx)
{
  int dimIt;
  int ptIt;

  int ndims = simpx->ndims;
  int npts  = simpx->npts;

  /*** Print Header ************************************/
  printf ("==========================================\n");
  printf ("  i");
  for (dimIt = 0; dimIt < ndims; dimIt++)
    {
      printf (" points[%d]", dimIt);
    }
  printf (" vals[i]\n\n");

  /************ Print Data ************************************/
  for (ptIt = 0; ptIt < npts; ptIt++)
    {
      printf ("%3d ", ptIt);
      for (dimIt = 0; dimIt < ndims; dimIt++)
	{
	  printf ("%9.4f ", simpx->points[ptIt * ndims + dimIt]);
	}
      printf ("%9.6f\n", simpx->vals[ptIt]);
    }
  printf("\n\n");

  return;
}


/***************************************************************************//**
 * check_borders() changes ptry to be within simpx's borders if necessary
 *
 * @param simpx 
 *     Address of the simplex whose borders we want to make sure ptry respects.
 *     Neither the address simpx points to nor the simplex simpx points to 
 *     change.
 * @param ptry 
 *     Address of the point to be checked.  The function may change one or
 *     more components of ptry to ensure the boundaries of simpx are respected.
 *  
 * @return number of times the borders of simpx are crossed
*******************************************************************************/
int check_borders( struct simplex const * const simpx, double * const ptry)
{
  int crossings = 0;
  int dimIt;
  int const ndims = simpx->ndims;
  double noise;
  double const * const borders = simpx->borders;

  if( borders == NULL )
  {
    return 0;
  }

  for( dimIt = 0; dimIt < ndims; dimIt++ )
  {
    noise = (borders[ 2*ndims + dimIt ] - borders[ 1*ndims + dimIt ]) 
            * .0010 * fabs( sin( 2.20 * crossings ));
    assert( noise >= 0 );
    /* Check crossing of hard lower borders */
    if ( ( simpx->borders[0*ndims+dimIt] == 0 || 
	   simpx->borders[0*ndims+dimIt] == 2 ) && 
	   ptry[dimIt] < simpx->borders[1*ndims+dimIt] )
    {
      ptry[dimIt] = borders[1*ndims + dimIt] + noise;
      crossings++;
    }

    /* Check crossing of hard upper borders */
    if ( ( simpx->borders[0*ndims+dimIt] == 0 || 
	   simpx->borders[0*ndims+dimIt] == 1 ) && 
	 ptry[dimIt] < simpx->borders[2*ndims+dimIt] )
    {
      ptry[dimIt] = borders[2*ndims + dimIt] - noise;
      crossings++;
    }

  }
      
  return crossings;

}


/***************************************************************************//**
 * Replaces the highest point in simpx with ptry and updates psum & vals[ihi]
 *
 * @param simpx 
 *     Address of simplex to be updated
 * @param ptry
 *     point to replace simpx->points[simpx->ihi] with.  The function does
 *     not check that ptry is within the hard boundaries.
 *
*******************************************************************************/
void replace_pt (
    struct simplex * const simpx, 
    double const * const ptry )
{
  int ndims   =  simpx->ndims;
  int ihi     =  simpx->ihi;
  int dimIt;
  double ytry = ptry[ndims];

  /* Replace highest point with ptry, save change to point, and update psum */
  simpx->vals[ihi] = ytry;
  for( dimIt = 0; dimIt < ndims; dimIt++)
  {
    simpx->psum[dimIt] += ptry[dimIt] - simpx->points[ihi * ndims + dimIt];
    simpx->change[dimIt] = simpx->points[ihi*ndims + dimIt] - ptry[dimIt];
    simpx->points[ihi * ndims + dimIt] = ptry[dimIt];
  }

  return;
}
