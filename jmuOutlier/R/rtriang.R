rtriang <-
function(n, min=0, max=1) {
  # Random generation for the symmetric triangular distribution with endpoints \code{min} and \code{max}.
  # 'n': Number of observations. If code{length(n)>1}, the length is taken to be the number required.
  # 'min': Left endpoint of the triangular distribution.
  # 'max': Right endpoint of the triangular distribution.
  # example:   # 20 random variates from a Triangular( 100, 200 ) distribution.
  #            rtriang( 20, 100, 200 )
  if (!is.numeric(n))  stop("'n' must be numeric.")
  if (length(n)>1)  n=length(n)
  if (n<1)  stop("'n' must be a positive integer.")
  if (!is.numeric(min))  stop("'min' must be numeric.")
  if (!is.numeric(max))  stop("'max' must be numeric.")
  if (min==max) stop("Endpoints cannot be equal.")
  if (min>max) {temp=min; min=max; max=temp};    ( runif(n,min,max) + runif(n,min,max) ) /2
}
