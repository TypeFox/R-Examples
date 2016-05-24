ptriang <-
function(q, min=0, max=1, lower.tail=TRUE) {
  # Cumulative distribution function for the symmetric triangular distribution with endpoints \code{min} and \code{max}.
  # 'q': Vector of quantiles.
  # 'min': Left endpoint of the triangular distribution.
  # 'max': Right endpoint of the triangular distribution.
  # 'lower.tail': Logical; if \code{TRUE} (default), probabilities are \code{P[X <= x]}; otherwise, \code{P[X > x]}.
  # example:    ptriang( seq( 100, 200, length.out=11 ), 100, 200 )
  # example:    ptriang( seq( 100, 200, length.out=11 ), 100, 200, FALSE )
  if (!is.numeric(q))    stop("'q' must be numeric.")
  if (!is.numeric(min))  stop("'min' must be numeric.")
  if (!is.numeric(max))  stop("'max' must be numeric.")
  if (!is.logical(lower.tail))  stop("'lower.tail' must be logical.")
  if (min==max) stop("Endpoints cannot be equal.")
  if (min>max) {temp=min; min=max; max=temp} ;    c=(min+max)/2
  p=(q-min)^2/(max-min)/(c-min)*(min<=q&q<=c)+(1-(max-q)^2/(max-min)/(max-c))*(c<q&q<=max)+(q>max)
  if (!lower.tail) p=1-p;    return(p)
}
