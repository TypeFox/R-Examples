rlaplace <-
function(n, mean=0, sd=1) {
  # Laplace (double exponential) random generation with mean equal to \code{mean} and standard deviation equal to \code{sd}. 
  # 'n': Number of observations. If code{length(n)>1}, the length is taken to be the number required.
  # 'mean': Population mean.
  # 'sd': Population standard deviation.
  # example:   # 20 random variates from a Laplace( 50, 10 ) distribution.
  #            rlaplace( 20, 50, 10 )
  if (!is.numeric(n))  stop("'n' must be numeric.")
  if (length(n)>1)  n=length(n)
  if (n<1)  stop("'n' must be a positive integer.")
  if (!is.numeric(mean))  stop("'mean' must be numeric.")
  if (!is.numeric(sd))  stop("'sd' must be numeric.")
  if (sd<0)  stop("'sd' cannot be negative.")
  rexp(n, (1/sd)) * sample(c(-1, 1), n, replace=TRUE) / sqrt(2) + mean
}
