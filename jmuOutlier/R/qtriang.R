qtriang <-
function(p, min=0, max=1) {
  # Quantile function for the symmetric triangular distribution with endpoints \code{min} and \code{max}.
  # 'p': Vector of probabilities.
  # 'min': Left endpoint of the triangular distribution.
  # 'max': Right endpoint of the triangular distribution.
  # example:  # 5th, 15th, 25th, ..., 95th percentiles from a Triangular( 100, 200 ) distribution.
  #           qtriang( seq( 0.05, 0.95, length.out=11 ), 100, 200 )
  if (!is.numeric(p))  stop("'p' must be numeric.")
  if (p<0 | p>1)  stop("'p' must be between 0 and 1.")
  if (!is.numeric(min))  stop("'min' must be numeric.")
  if (!is.numeric(max))  stop("'max' must be numeric.")
  if (min==max) stop("Endpoints cannot be equal.")
  if (min(p)<0 || max(p)>1) stop("Values of `p' must be between zero and one.")
  if (min>max) {temp=min; min=max; max=temp}
  min + (max-min)*(  (0<=p&p<=0.5)*(0.5*sqrt(2*p))+(0.5<p&p<=1)*(0.5+0.5*sqrt(2*abs(p-0.5)))  )
}
