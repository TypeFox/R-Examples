plaplace <-
function(q, mean=0, sd=1, lower.tail=TRUE) {
  # Laplace (double exponential) cumulative distribution function with mean equal to \code{mean} and standard deviation equal to \code{sd}. 
  # 'q': Vector of quantiles.
  # 'mean': Population mean.
  # 'sd': Population standard deviation.
  # 'lower.tail': Logical; if \code{TRUE} (default), probabilities are \code{P[X <= x]}; otherwise, \code{P[X > x]}.
  # example:    plaplace( seq( 20, 80, length.out=11 ), 50, 10 )
  # example:    plaplace( seq( 20, 80, length.out=11 ), 50, 10, FALSE )
  if (!is.numeric(q))  stop("'q' must be numeric.")
  if (!is.numeric(mean))  stop("'mean' must be numeric.")
  if (!is.numeric(sd))  stop("'sd' must be numeric.")
  if (sd<0)  stop("'sd' cannot be negative.")
  if (!is.logical(lower.tail))  stop("'lower.tail' must be logical.")
  p=(q>mean) - ((q>mean)*2-1) * exp(-abs(q-mean)*sqrt(2)/sd)/2 
  if (!lower.tail) p=1-p;    return(p)
}
