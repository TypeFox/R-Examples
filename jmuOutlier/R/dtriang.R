dtriang <-
function(x, min=0, max=1) {
  # Probability density function for the symmetric triangular distribution with endpoints \code{min} and \code{max}.
  # 'x': Vector of quantiles.
  # 'min': Left endpoint of the triangular distribution.
  # 'max': Right endpoint of the triangular distribution.
  # example:    dtriang( seq( 100, 200, length.out=11 ), 100, 200 )
  if (!is.numeric(x))    stop("'x' must be numeric.")
  if (!is.numeric(min))  stop("'min' must be numeric.")
  if (!is.numeric(max))  stop("'max' must be numeric.")
  if (min==max) stop("Endpoints cannot be equal.")
  if (min>max) {temp=min; min=max; max=temp} ;    c=(min+max)/2
  2*(x-min)/(max-min)/(c-min)*(min<=x & x<=c) + 2*(max-x)/(max-min)/(max-c)*(c<x & x<=max)
}
