#include <omp.h>
#include "simplex.h"

/***************************************************************************//**
 * Uses a variation on the Nelder-Mead simplex algorithm from Numerical Recipes
 *
 * This function uses the same logic as the code in section 10.5 of Numerical
 * Recipes but was rewritten so be more clear and to be parallelizable.
 *
 * The four simplex actions are each run in their own thread:
         - reflect (the highest point across the opposite face)
         - extend ( /reflec
         - contract
         - reflect+contract
 * The threads then rejoin and one thread is used to determine which new node
 * should replace the simplex's high point.  If none of the 4 possibilities 
 * are improvments then the main thread performs an all contraction serially.
 *
 * FIXME -- Crossings code needs tested
 * FIXME -- 3*NVD Magic Number: Should this be an argument?  Why 3*?
 * FIXME -- Should NMAX be an argument?
 *
 * param simpx 
 *     The simplex which will be used to find a minimum of the function funk
 * param ftol
 *     Function returns when the ratio of the range of function values to the 
 *     sum of the highest and lowest values is less than half this value
 * param funk
 *     *funk() is the cost function that each point of the simplex will be
 *     evaluated at.  Its arguments are the dimension of the simplex and
 *     a pointer to the paramaters at which it should be evaluated.
 *
 * return the index of the lowest point in the simplex function
*******************************************************************************/ 
int amoeba_omp (
    struct simplex * const simpx, 
    double const ftol, 
    double (*funk)(int, double const * const), 
    int * const nfunk,
    double const momentum
    )
{
  double *vals = simpx->vals;
  int *ihi = &simpx->ihi;
  int *inhi = &simpx->inhi;
  int *ilo = &simpx->ilo;
  int ptIt, dimIt;
  double rtol;

  int totalCrossings = 0;

  double *Reflect     = malloc ( sizeof (double) * (simpx->ndims+1));
  int    xReflect;	//Number of times the Reflection point crossed borders
  double *Extend      = malloc ( sizeof (double) * (simpx->ndims+1));
  int    xExtend;
  double *RefContract = malloc ( sizeof (double) * (simpx->ndims+1));
  int    xRefContract;
  double *Contract    = malloc ( sizeof (double) * (simpx->ndims+1));
  int    xContract;
  assert ( Reflect && Extend && RefContract && Contract );
  
  get_psum( simpx );  /* FIXME call may be better in initialize_simplex() */

  for( *nfunk=0;  *nfunk  < NMAX; (*nfunk)+=4)
  {

    find_hilo(simpx);

    /**** Check for convergance or excessive crossings ****/
    rtol = 2.0 * fabs(      vals[*ihi]  -      vals[*ilo] )/ 
                     ( fabs(vals[*ihi]) + fabs(vals[*ilo]) + 1e-323);
                     /* Adding 1e-323 to ensure divisor > 0 */

    if (rtol < ftol || totalCrossings > 3*NVD )
    {
      free(Reflect);
      free(Extend);
      free(RefContract);
      free(Contract);
      return *ilo;
    } 

#ifdef DEBUG
    print_simplex(simpx);
#endif


    /**** Try to improve the worst point in simplex ****/
    /* Calculate potential points and function value */
    /* This single #pragma directive causes 4 threads to be created as 
     * specified by the argument to num_threads().  The arguments shared() 
     * get passed to each thread by reference.  Care is necessary to prevent
     * multiple threads from writing to the same variable at the same time.
     * Each thread then executes the following code block concurrently. */
#pragma omp parallel num_threads(4) default(none)\
    shared (funk,  Reflect,  Extend,  RefContract,  Contract, \
	          xReflect, xExtend, xRefContract, xContract )
    {
      int thread_id = omp_get_thread_num(); 
      switch(thread_id)
      {
	case 0:
	  xReflect     = amotry( simpx, funk, -1.0, Reflect,     momentum);
	  break;
	case 1:
	  xExtend      = amotry( simpx, funk, -2.0, Extend,      momentum);
	  break;
	case 2:
	  xRefContract = amotry( simpx, funk, -0.5, RefContract, momentum);
	  break;
	case 3:
	  xContract    = amotry( simpx, funk,  0.5, Contract,    momentum);
	  break;
      }
    } /* Threads return control to single process here */

#ifdef DEBUG
    printf("Reflect: %lf; Extend: %lf; RefCont: %lf; Contract %lf\n",
	Reflect[simpx->ndims], Extend[simpx->ndims], RefContract[simpx->ndims], Contract[simpx->ndims]);
#endif
    
    /* Reflect or Extend */
    if( Reflect[simpx->ndims] < simpx->vals[*ilo] )
    {
      if ( Extend[simpx->ndims] < Reflect[simpx->ndims] )
      {
	replace_pt( simpx, Extend );
	totalCrossings += xExtend;
#ifdef DEBUG
	puts("Extended");
#endif
      }
      else
      {
	replace_pt( simpx, Reflect );
	totalCrossings += xReflect;
#ifdef DEBUG
	puts("Reflected");
#endif
      }
    }
    /* Reflect */
    else if( Reflect[simpx->ndims] < simpx->vals[*inhi] )
    {
	replace_pt( simpx, Reflect );
	totalCrossings += xReflect;
#ifdef DEBUG
	puts("Reflected");
#endif
    }
    /* Reflect or Reflect&Contract */
    else if( Reflect[simpx->ndims] < simpx->vals[*ihi] )
    {
      if ( RefContract[simpx->ndims] < Reflect[simpx->ndims] )
      {
	replace_pt( simpx, RefContract );
	totalCrossings += xRefContract;
#ifdef DEBUG
	puts("RefContracted");
#endif
      }
      else
      {
	replace_pt( simpx, Reflect );
	totalCrossings += xReflect;
#ifdef DEBUG
	puts("Reflected");
#endif
      }
    }
    /* Contract or All Contract */
    else /* Reflect[simpx->ndims] >= simpx->vals[*ihi] */
    {
      if ( Contract[simpx->ndims] < simpx->vals[*ihi] )
      {
	replace_pt( simpx, Contract );
	totalCrossings += xContract;
#ifdef DEBUG
	puts("Contracted");
#endif
      }
      else /* Contract All towards best point */
      {
#pragma omp parallel for private(dimIt)
	for( ptIt = 0; ptIt < simpx->npts; ptIt++ )
	{
	  if ( ptIt == *ilo )
	    continue;

	  for( dimIt = 0; dimIt < simpx->ndims; dimIt++ )
	  {
	    simpx->points[ptIt * simpx->ndims + dimIt] = 
	      0.5 * (simpx->points[ptIt * simpx->ndims + dimIt] + 
		     simpx->points[*ilo * simpx->ndims + dimIt]);
	  }
	  vals[ptIt] = (*funk)(simpx->ndims, &(simpx->points[ptIt*simpx->ndims]));
	}

	*nfunk += simpx->ndims;
	get_psum(simpx);
#ifdef DEBUG
	puts("All Contracted");
#endif
      }
    }
	  
  }

  /*** Too many function evals ***/
  free(Reflect);
  free(Extend);
  free(RefContract);
  free(Contract);
  return (*ilo);

}


/***************************************************************************//**
 * Creates a new point along the line from ihi through the opposite face
 *
 * Function averages all the points of the simplex except ihi to find the center
 * of the face opposite ihi.  It then finds the point along the line from ihi
 * through this centroid that corresponds to fac.  We bring the new point back
 * into the specified borders if it has crossed, calculate the function at this
 * point and return.
 * 
 *
 * @param simpx 
 *     The simplex which will be used to find a minimum of the function funk
 * @param funk
 *     *funk() is the cost function that each point of the simplex will be
 *     evaluated at.  Its arguments are the dimension of the simplex and
 *     a pointer to the paramaters at which it should be evaluated.
 * @param fac
 *     Determines the distance from the face opposite ihi of the new point.
 *     A point with fac=-1 will be the same distance from the face as ihi only
 *     in the opposite direction.
 * @param ptry
 *     pointer to simpx->ndims+1 doubles.  When the function returns the first 
 *     ndims doubles will contain the extrapolated point.  The last double
 *     will contain the function value at that point.
 *
 * @return number of borders ptry crossed
 *     
*******************************************************************************/ 
int amotry (
    struct simplex const * const simpx, 
    double (*const funk)(int const, double const * const), 
    double const fac, 
    double * const ptry,
    double const momentum )
{
  int dimIt;
  double p_bar;
  int crossings;

  /* Create aliases for frequently used elements of simpx */
  int ndims   =  simpx->ndims;
  int ihi     =  simpx->ihi;
  double *vals =  simpx->vals;

  assert( ndims );
  for (dimIt = 0; dimIt < ndims; dimIt++)
  {
    /* Find centroid of the points < vals[ihi] in the current dimension */
    p_bar = ( simpx->psum[dimIt] - simpx->points[ ihi * ndims + dimIt ] ) 
            / ndims;

    /* Reflect points[ihi] about p_bar by a factor of "fac" */
    ptry[dimIt] = ( simpx->points[ihi * ndims + dimIt] * fac ) + 
                  ((1-fac)*p_bar) + momentum*simpx->change[dimIt] ;
  }

  crossings = check_borders( simpx, ptry ); 
  
  
  ptry[simpx->ndims] = (*funk)(simpx->ndims, ptry);

  return crossings;  /* Return number of crossings */

}

