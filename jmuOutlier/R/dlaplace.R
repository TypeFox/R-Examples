dlaplace <-
function(x, mean=0, sd=1) {
  # Laplace (double exponential) density function with mean equal to \code{mean} and standard deviation equal to \code{sd}. 
  # 'x': Vector of quantiles.
  # 'mean': Population mean.
  # 'sd': Population standard deviation.
  # example:    dlaplace( seq( 20, 80, length.out=11 ), 50, 10 )
  if (!is.numeric(x))  stop("'x' must be numeric.")
  if (!is.numeric(mean))  stop("'mean' must be numeric.")
  if (!is.numeric(sd))  stop("'sd' must be numeric.")
  if (sd<0)  stop("'sd' cannot be negative.")
  exp(-abs(x-mean)*sqrt(2)/sd)/(sd*sqrt(2))
}
