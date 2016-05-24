#' @title Expectation of CUB models
#' @description Compute the expectation of a CUB model without covariates.
#' @aliases expcub00
#' @usage expcub00(m, pai, csi)
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter
#' @export expcub00
#' @seealso \code{\link{varcub00}}, \code{\link{expcube}}, \code{\link{varcube}}
#' @keywords distribution
#' @references Piccolo D. (2003). On the moments of a mixture of uniform and shifted binomial random variables. 
#' \emph{Quaderni di Statistica}, \bold{5}, 85--104
#' @examples
#' m<-10
#' pai<-0.3
#' csi<-0.7
#' meancub<-expcub00(m,pai,csi)


expcub00 <-
function(m,pai,csi){(m-1)*pai*(0.5-csi)+(m+1)/2}
