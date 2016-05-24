qlaplace <-
function(p, mean=0, sd=1) {
  # Laplace (double exponential) quantile function with mean equal to \code{mean} and standard deviation equal to \code{sd}. 
  # 'p': Vector of probabilities.
  # 'mean': Population mean.
  # 'sd': Population standard deviation.
  # example:  # 5th, 15th, 25th, ..., 95th percentiles from a Laplace( 50, 10 ) distribution.
  #           qlaplace( seq( 0.05, 0.95, length.out=11 ), 50, 10 )
  if (!is.numeric(p))  stop("'p' must be numeric.")
  if (p<0 | p>1)  stop("'p' must be between 0 and 1.")
  if (!is.numeric(mean))  stop("'mean' must be numeric.")
  if (!is.numeric(sd))  stop("'sd' must be numeric.")
  if (sd<0)  stop("'sd' cannot be negative.")
  mean + ( 2*(p<=0.5)-1 ) * sd / sqrt(2) * log( 2 * (p*(p<=0.5)+(1-p)*(p>0.5)) )
}
