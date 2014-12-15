#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include "simplex.h"

double gk_WClockSeconds(void);
#define gk_clearwctimer(tmr) (tmr = 0.0 )
#define gk_startwctimer(tmr) (tmr -= gk_WClockSeconds() )
#define gk_stopwctimer(tmr)  (tmr += gk_WClockSeconds() )
#define gk_getwctimer(tmr)   (tmr)

/* A 2D Polynomial function with minimum of 10 at (-3,0)  *******************/
double poly_2d ( int const, double const * const );
/* A 2D Polynomial function with minimum of 0 at (-3,0)  *******************/
double poly_2d_0 ( int const, double const * const );

/***** Main Driver Function ****************************/
int main (void)
{
  struct simplex *simpx;
  int best;
  int nfunc; 
  double answer[5] = {-3.0, 0.0};

  int ndims = 2;
  int npts  = 3;
  int i, success = 0;

  /** Check 2D Polynomial ***************************************/
  double *pts = malloc ( sizeof(double) * ndims*npts ); /* Freed in destroy_simplex */
  /* Set initial Points */
  pts[0] = -3.9; pts[1] =  0.1; 
  pts[2] =  3.5; pts[3] =  1.0;
  pts[4] =  5.0; pts[5] = -9.0;
  
  /* Set Soft Boundaries */
  double *bounds = malloc ( sizeof(double) * 3*ndims ); /* Freed in destroy_simplex */
  bounds[0] =   3; bounds[1] =   3;
  bounds[2] = -10; bounds[3] = -10;
  bounds[4] =  10; bounds[5] =  10;

  nfunc = 0;
  simpx = initialize_simplex(2, 3, 2, pts, bounds, poly_2d);
  assert( simpx );
  best = amoeba_omp (simpx, 1.0e-6, poly_2d, &nfunc, 0);
  printf("\n\n\n");
  printf("Poly_2d: Best Point: ");
  for ( i = 0; i < simpx->ndims; i++)
    {
      printf("%.3f  ", simpx->points[best*simpx->ndims + i]);
      if( fabs( answer[i] - simpx->points[best*simpx->ndims + i] ) < .01 )
        {
          success++;
        }
    }
  printf("\n\tMinimum at Best Point: %.3f ", simpx->vals[best]);
  printf("\n\tFunction Calls: %d\n\t", nfunc);
  success == simpx->ndims ? printf("SUCCEEDED\n\n") : printf("FAILED\n\n");
  destroy_simplex( simpx );
  /********************************************************************/


  /** Check 2D Polynomial with minimum of 0.   ***************************/ 
  /** Verifies no divide by 0 error.   ***********************************/
  nfunc = 0;
  double *bound = malloc ( sizeof(double) * 3*ndims ); /* Freed in destroy_simplex */
  bound[0] =   3; bound[1] =   3;
  bound[2] = -10; bound[3] = -10;
  bound[4] =  10; bound[5] =  10;
  simpx = initialize_simplex(2, 3, 1, NULL, bound, poly_2d_0);
  best = amoeba_omp (simpx, 1.0e-6, poly_2d_0, &nfunc, 0);
  success = 0;



  printf("Poly_2d_0: Best Point: ");
  for ( i = 0; i < simpx->ndims; i++)
    {
      printf("%.3f  ", simpx->points[best*simpx->ndims + i]);
      if( fabs( answer[i] - simpx->points[best*simpx->ndims + i] ) < .01 )
        {
          success++;
        }
    }
  printf("\n\tMinimum at Best Point: %.3f ", simpx->vals[best]);
  printf("\n\tFunction Calls: %d\n\t", nfunc);
  success == simpx->ndims ? printf("SUCCEEDED\n\n") : printf("FAILED\n\n");
  destroy_simplex( simpx );
  /********************************************************************/

  return 0;
}

double gk_WClockSeconds(void)
{
  struct timeval ctime;

  gettimeofday(&ctime, NULL);

  return (double) ctime.tv_sec + (double).000001*ctime.tv_usec;
}

/* A 2D Polynomial function with minimum of 10 at (-3,0)  *******************/
double poly_2d (int const n, double const * const x)
{
  return 3*(x[0]+3)*(x[0]+3)  +  4*(x[1]-0)*(x[1]-0)  +  10;
}

/* A 2D Polynomial function with minimum of 0 at (-3,0)  *******************/
double poly_2d_0 (int const n, double const * const x)
{
  return 3*(x[0]+3)*(x[0]+3)  +  4*(x[1]-0)*(x[1]-0)  +  0;
}
