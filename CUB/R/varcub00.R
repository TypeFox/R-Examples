#' @title Variance of CUB models without covariates
#' @description Compute the variance of a CUB model without covariates.
#' @aliases varcub00
#' @usage varcub00(m, pai, csi)
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter 
#' @export varcub00
#' @seealso  \code{\link{CUB}}, \code{\link{expcub00}}, \code{\link{probcub00}}
#' @references 
#' Piccolo D. (2003). On the moments of a mixture of uniform and shifted binomial random variables. 
#' \emph{Quaderni di Statistica}, \bold{5}, 85--104
#' @keywords distribution
#' @examples 
#' m<-9
#' pai<-0.6
#' csi<-0.5
#' varcub<-varcub00(m,pai,csi)

varcub00 <-function(m,pai,csi){
  (m-1)*(pai*csi*(1-csi)+(1-pai)*(((m+1)/12) + pai*(m-1)*(csi-1/2)^2  ))
}
