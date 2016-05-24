##' Function to calculate statistical moments
##' 
##' The function calculates the first 4 moments, i.e. the mean, variance, skew,
##' kurtosis.
##' 
##' The units of the first moment are the same as x, the units of the second
##' moment are x\eqn{\mbox{\textasciicircum}}{^}2, and the third and fourth
##' moments are dimensionless.
##' 
##' @param count A vector of the observed instances per class
##' @param x A vector of the same length as count defining the class. If
##' missing, and if count is of class spectral, then x is equal to
##' trackfreq(count). If x is missing and is not of class spectral, then x
##' default to 0:(length(count)-1)
##' @param minval If T, subtract min(count) from count so that the minimum
##' value of count is zero. This is principally used in calculating spectral
##' moments where count is in decibels, and more generally if count contains
##' negative values.
##' @author Jonathan Harrington
##' @references Snedecor, G & Cochran, W. 'Statistical Methods' Iowa State
##' Press. Wuensch,K., 2005
##' @keywords math
##' @examples
##' 
##' # first four moments of a vector
##' mom <- moments(bridge[,2])
##' # the above is the same as moments(bridge[,2], 0:12)
##' # first four moments of a spectral vector with the dB values
##' # reset so that the minimum dB value is 0. The d.c. offset is also
##' # excluded in the calculation
##' mom <- moments(e.dft[-1], minval=TRUE)
##' # the temporal skew of F1 for the 10th segment. Use
##' m <- moments(vowlax.fdat[10,1]$data)[3]
##' 
##' 
##' @export moments
"moments" <- function(count, x, minval = FALSE)
  
{
  # compute moments. x is a numeric class
  # count is the frequency with which that 
  # particular class occurs
  # This function gives exactly the same
  # results as those for the mean, variance
  # skewness and kurtosis in example Table 3.13.1
  # p. 87, Snedecor & Cochran, 'Statistical Methods'
  # 6th Edition, 1975. Let the arguments count and x
  # equal f and U respectively in their example
  # the centre of gravity with minval = F.
  # the first two moments in this function
  # also give the same results as in Harrington & Cassidy.
  if(minval)
    count <- count - min(count)
  if(missing(x))
  {
    if(is.spectral(count))
      x <- trackfreq(count)
    else
      x <- 0:(length(count)-1)
  }
  k <- 1
  mom1 <- sum((x - 0)^k * count) / sum(count)
  # the variance
  k <- 2
  mom2 <- sum((x - mom1)^k * count) / sum(count)
  
  # third moment
  k <- 3
  mom3 <- (sum((x - mom1)^k * count) / sum(count)) / (mom2 * sqrt(mom2))
  
  # fourth moment
  k <- 4
  # peaked distributions show positive kurtosis
  # flat-topped distributions show negative kurtosis
  mom4 <- (sum((x - mom1)^k * count) / sum(count)) / mom2^2 - 3
  c(mom1, mom2, mom3, mom4)
}
