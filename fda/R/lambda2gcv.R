lambda2gcv <- function(log10lambda, argvals, y, fdParobj, 
                          wtvec=rep(1,length(argvals)))
{
#  LAMBDA2GCV smooths data using smooth_basis using 
#  LAMBDA = 10^LOG10LAMBDA, and returns the GCV value that result. 

#  Return:
#  GCV ... a vector or matrix of GCV values, depending on whether Y is a
#          matrix of array, respectively. 

#  Last modified 21 October 2008 by Jim Ramsay

  fdParobj$lambda <- 10^log10lambda
  gcv <- smooth.basis(argvals, y, fdParobj, wtvec)$gcv

  return(gcv)
}