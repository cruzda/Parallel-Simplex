#ifndef __DB_SIMPLEX_H__
#define __DB_SIMPLEX_H__

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define NMAX 100000
#define NVD  50

struct simplex { 
    int ndims;     /* Number of dimensions in problem */
    int npts;      /* Number of points in simplex.   */

/**points:
*    Points to a block of memory of size >= npts*ndims*sizeof(double) and 
*    contains the points of the simplex.
*    The memory block can be thought of as a two dimensional array where each 
*    "row" begins at the index [row*ndims]. Each row contains ndims doubles 
*    specifing the location of one point in the simplex.                    */
    double * points; 

/**change:
 *   Contains the magintude and direction of the last change to simplex using 
 *   replace_pt(); */
    double * change;
    
/**borders:
*    Points to a block of memory of size >= 3*ndims*sizeof(double) which contains 
*    information defining two boundaries for the simplex.  A soft boundary which
*    is used when source = 1 to determine where the initial simplex may be  
*    located.  The soft boundary range is also used when bringing the simplex 
*    back within a hard boundary.  The hard boundary is not to be crossed by  
*    the simplex at all.  Boundaries are not strictly enforced and are only 
*    checked when the function check_borders() is called.  If borders need to be 
*    enforced make sure to call check_borders() before calling replace_pt(). 
*
*    Boundaries are specified using three "rows" of ndims each.   
*    The first row which starts at borders[0] is expected to contain integer  
*    values of 0, 1, 2 or 3 (how the program handles other values is  
*    undetermined).  These values have the following meanings: 
*        0:  Both the high and low values are hard boundaries 
*        1:  Only the low value is a hard boundary.  The high value is soft. 
*        2:  Only the high value is a hard boundary.  The low value is soft. 
*        3:  Neither value is a hard boundary.  Both are soft. 
*    The second "row" begins at borders[ndims] and is expected to contain doubles
*    specifying the lower boundaries.   
*    The third "row" begins at borders[ndims*2] and is expected to contain  
*    doubles specifying the upper boundaries. These should be greater than the 
*    corresponding lower boundary but the program doesn't verify this.  
*    borders may also be NULL in which case no borders are enforced.*/
    double * borders; 

    /* Array containing value of function at corresponding point in points */
    double * vals; 

    /* (psum - points[ihi])/ndims = centroid of face opposing the 
     *                              highest point  */
    double * psum; 

 
    int ilo;      /* index to point with lowest function value */
    int inhi;     /* index to point with 2nd highest function value */
    int ihi;      /* index to point with highest function value */
}; 
 
/* Allocates memory for simplex and initializes */
struct simplex * initialize_simplex(
    int const ndims,
    int const npts, 
    int const source, 
    double * const pts,
    double * bounds,
    double (* const func)(int const, double const * const));

/* Finds the point obtained by adding each point in points as if it were a vector */
void get_psum(struct simplex * const simpx);

/* Finds the highest, 2nd highest, and lowest points in points array */
void find_hilo(struct simplex * const simpx);
int check_border(struct simplex const * const simpx, double * const ptry);

/* replaces point ihi in simplex with ptry updating psum and vals */
void replace_pt (
    struct simplex * const simpx, 
    double const * const ptry);

/* Prints points of the simplex and their function values */
void print_simplex(struct simplex const * const simpx); 

/* Frees memory allocated in initialize simplex(); */
void destroy_simplex(struct simplex * simpx); 

/* Creates a new point along the line from ihi through the oppose face */
int amotry (
    struct simplex const * const simpx, 
    double (* const funk) (int const, double const * const), 
    double const fac, 
    double * const ptry,
    double const momentum);

/* Uses a varion on Nelder-Mead to find a minimum of a function using 4 cores */
int amoeba_omp (
    struct simplex * const simpx, 
    double const ftol, 
    double (* const funk) (int const, double const * const), 
    int * const nfunk,
    double const momentum );

#endif
