BJnoint <- function(x, y, delta, beta0 = NA, maxiter = 30, error = 0.00001)
{
# This is an R function to compute the Buckley-James estimator for
# censored regression WITHOUT intercept term. (but you may still force x 
# to have a column of 1's). It calls another function iter(). This function
# use to call a C function.  Written by mai zhou, mai@ms.uky.edu
# First (C) version Jan. 1994. Last C revision: May 10, 1995
# R version: 2004 Jan 20. R speed actually is OK, no need of C function.
#
# Input:
# x is a matrix of N rows (covariates).
# y is the observed (censored) responses --- a vector of length N.
# delta is a vector of length N. delta =1 means (y) is not censored.
#           delta = 0 means y is right censored, i.e. the true
#        response is larger than y.
#
# Optinal input: maxiter = number of maximum iteration, 
#                           default is 30.  minimum iteration 
#                           is internally set to 3.
#                error = when the consecutive iteration changes less then
#                        error, the iteration will stop. default is .00001
#                beta0 = initial estimator.
# Output:
# the estimate, beta, and an extra integer at the end: number of iterations.
#
# Bug: More careful for ties. Do something better when last obs. is censored.

     x <- as.matrix(x)            # to insure x is a matrix, not a vector
     newtemp <- matrix(NA, ncol = ncol(x), nrow = 3)
     newtemp[1,] <- beta0
     if(any(is.na(beta0))) 
        newtemp[1, ] <- lm(y ~ 0 + x)$coef  # get initial est.
     for(i in 2:3) {
           newtemp[i, ] <- iter(x, y, delta, newtemp[i - 1, ])
     }
     num <- 2                        # do at least 2 iterations
     while(num <= maxiter && error <= sum(abs(newtemp[2, ] - newtemp[3, ])))
          {
          newtemp[2, ] <- newtemp[3, ] 
          newtemp[3, ] <- iter(x, y, delta, newtemp[2, ])
          num <- num + 1
          }
    if(num > maxiter)   # always take average? or only when not converge?
    {newtemp[3, ] <- (newtemp[2, ] + newtemp[3, ])/2} 
    list(beta=newtemp[3, ], iteration=num)
}
