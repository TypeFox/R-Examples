#' Generate independent Brownian bridges. If \code{BM} is specified
#' then bridges are constructed from a given Brownian Motion \code{BM}
#' \deqn{Y_t(u) = BM_t(u) - u BM_t(1)}
#' for each \eqn{t} and \eqn{u \in [0,1]}.
#' Otherwise \code{n} Brownian bridges are generated indipendently.
#'
#' @title Generate brownian bridges
#' @param n number of observations to generate
#' @param BM brownian motion to use
#' @param d dimension (sampling at d points)
#' @return n x d matrix with independent n observations
#' @export
#' @importFrom graphics plot
#' @examples
#' bm = rbm(100)
#' plot(bm)
#' @seealso \code{\link{rbm}}
rbb = function(n=NULL,d=100,BM=NULL)
{  
  if (is.null(BM) && is.null(n))
    stop("Either a Brownion motion BM or number of oservations n should be provided.")
  if (!is.null(n) && !is.positiveint(n))
    stop("n must be a positive integer.")
  
  if (is.null(BM)){
    BM = rbm(n,d)
  }
  else {
    d = dim(BM)[2]
  }
  Rt = BM[,d] %*% t(1:d/d)
  BM - Rt
}
