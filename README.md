# Parallel Simplex

This repository contains an implementation of Nelder-Mead's Downhill Simplex
Method.  The algorithm has been modified to allow the use of 4 cores when 
calculating the various possible changes to the simplex at each step.

Optimizing a function involves two steps.  The first is to create a simplex
by calling the initialize_simplex() function.  It's arguments are described
in simplex.h and contain:
- the simplex's initial points (optional)
- boundaries for the simplex (optional)
- the cost function that needs to be optimized.

Pass newly created simplex to amoeba_omp() which will return the best point
found by the algorithm.

The simplex's current points can be printed with print_simplex( ).

Boundaries have not been throughly tested.  Use with caution.

To prevent memory leaks, be sure to destroy the created simplex with 
destroy_simplex().



##Testing

Code can be tested against two polynomial functions: 

3(a+3)^2 + 4(b-0)^2 + 10 ==> min of 10 @(-3,0)
3(a+3)^2 + 4(b-0)^2 +  0 ==> min of  0 @(-3,0)

by running 

$ make test && ./test



This simplex code was tested against the functions at
http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml starting with
the easiest.

To create an initial simplex I provided the first point from the data file
then created the remaining n points by scaling one of the n dimensions
by a constant factor.

I did not find any problems with the code itself.  Occasionally an error
will rear its head but it was never something that is a problem with
amoeba or that amoeba could detect.

The algorithm successfully found correct minimums for each of the first
5 functions (misra1a, chwirut2, chwirut1, Lanczos3 and Gauss1).

- On Gauss2, the algorithm failed. It found minimums 2-3x larger than
the certified minimum when starting with a simplex based on the provided
starting values.  Starting with the certified arguments resulted in
finding the correct minimum.  I also tried to use momentum to improve
the results but continued to find the same minimum regardless of the
size of the momentum (all 100 values between 0.00 & 1.00)

-   Algorithm was successful with Kirby2 w/o momentum.

-   It failed to find the same minimum for Hahn1, though the minimum
it found was within 1 standard deviation of the certified value.  I did
not run this with momentum.

-   On MGH17 it found minimums of 1.023 and larger using momentum up to
~.35, but the certified minimum was 5.5e-5.

-   In the Lanczos1 algorithm it did not find the certified minimum,
though momentum of .14 was helpful in reducing the minimum found by at
least 3 orders of magnitude.

-   Without momentum Lanczos2 did not find the certified minimum, but
with momentum = .03 &.04 among other values it did, though the simplex
didn't converge and amoeba returned after calling func 100000x.
