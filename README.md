# Parallel Simplex

##

I've been creating code to test the simplex against the functions at http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml starting with the easiest.

To create an initial simplex I provided the first point from the data file then created the remaining n points by scaling one of the n dimensions by a constant factor.

I have yet to find any problems with the code itself.  Occasionally an error will rear its head but it has yet to be something that is a problem with amoeba or that amoeba could detect.

The algorithm successfully found correct minimums for each of the first 5 functions (misra1a, chwirut2, chwirut1, Lanczos3 and Gauss1). 

On Gauss2 the algorithm failed finding minimums 2-3x larger than the certified minimum when starting with a simplex based on the provided starting values.  Starting with the certified arguments resulted in finding the correct minimum.  I also tried to use momentum to improve the results but continued to find the same minimum regardless of the size of the momentum (all 100 values between 0.00 & 1.00)

   Algorithm was successful with Kirby2 w/o momentum.

   It failed to find the same minimum for Hahn1, though the minimum it found was within 1 standard deviation of the certified value.  I did not run this with momentum.

   On MGH17 it found minimums of 1.023 and larger using momentum up to ~.35, but the certified minimum was 5.5e-5. 

   In the Lanczos1 algorithm it did not find the certified minimum, though momentum of .14 was helpful in reducing the minimum found by at least 3 orders of magnitude.

   Without momentum Lanczos2 did not find the certified minimum, but with momentum = .03 &.04 among other values it did, though the simplex didn't converge and amoeba returned after calling func 100000x.

   This was the extent of the functions I tested.  Attached is the code I've generated.  Please get in touch if you have any questions.
