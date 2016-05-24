# TODO: p-values for factors, etc.


#' Summarizing massively parallel linear model fits
#' 
#' \code{summary} method for class "\code{lm.mp}".
#' 
#' @method summary lm.mp
#' @param object an object of class \code{\link{lm.mp}}, ordinarily created by
#' the function of that name or by \code{\link{lm4d}}.
#' @param \dots not currently used.
#' @return \item{tstat}{matrix of pointwise t-statistics for each coefficient
#' in the linear model} \item{pvalue}{matrix of the pointwise p-values for each
#' coefficient in the linear model} \item{aicc}{vector of pointwise corrected
#' AIC}
#' @author Philip Reiss \email{phil.reiss@@nyumc.org} and Lei Huang
#' \email{huangracer@@gmail.com}
#' @seealso \code{\link{lm.mp}}
#' @examples
#' 
#' Y = matrix(rnorm(6000), nrow=20)
#' X = rnorm(20)
#' t1 = lm.mp(Y, ~X)
#' st1 = summary(t1)
#' @export
summary.lm.mp <- function(object, ...) {
	X = object$X
	coef = object$coef
	se.coef = object$se.coef
    n = dim(X)[1]
    p = dim(X)[2]
    tstat = coef[-1,] / se.coef[-1,]  
    pvalue = 2 * pt(-abs(tstat), n-p)  
	sigma2.mle = object$sigma2 * (n-p) / n
    aicc = n*log(sigma2.mle) + n*(n+p)/(n-p-2)	    
	list(tstat=tstat, pvalue=pvalue, aicc=aicc)
}

