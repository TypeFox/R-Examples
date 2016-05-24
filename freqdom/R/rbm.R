#' Generate \eqn{n} independent Brownian motions
#'
#' @title Generate brownian motions
#' @param n number of observations to generate
#' @param d dimension (sampling at d points)
#' @return n x d matrix with independent n observations
#' @importFrom graphics plot
#' @importFrom stats rnorm
#' @export
#' @examples
#' bm = rbm(100)
#' plot(bm)
#' @seealso \code{\link{rbb}}
rbm = function(n,d=100)
{  
  if (!is.positiveint(n))
    stop("n must be a positive integer.")

  sd = sqrt(1/(d-1))
  incr = matrix(rnorm((d-1)*n, mean=0, sd=sd),n,d-1)
  incr = cbind(rep(0,n),incr)
  t(apply(incr,1,cumsum))
}
