# Y: n x V outcome matrix (V=number of voxels, connections, etc.)
# formula: object of form "~ x1 + x2" (quotation marks not needed)


#' Massively parallel linear regression models
#' 
#' Efficiently fits \eqn{V} linear models with a common design matrix, where
#' \eqn{V} may be very large, e.g., the number of voxels in a brain imaging
#' application.
#' 
#' 
#' @param Y \eqn{n \times V} outcome matrix.
#' @param formula a formula object such as "\code{~ x1 + x2}".
#' @param store.fitted logical: Should the fitted values be stored?  For large
#' \eqn{V}, setting this to \code{TRUE} may cause memory problems.
#' @return \item{coef}{\eqn{p \times V} matrix of coefficient estimates.}
#' \item{sigma2}{\eqn{V}-dimensional vector of error variance estimates.}
#' \item{se.coef}{\eqn{p \times V} matrix of coefficient standard error
#' estimates.} \item{X}{\eqn{n \times p} common design matrix.}
#' \item{fitted}{\eqn{n \times V} matrix of fitted values.}
#' @author Philip Reiss \email{phil.reiss@@nyumc.org}, Lei Huang
#' \email{huangracer@@gmail.com}, and Yin-Hsiu Chen
#' \email{enjoychen0701@@gmail.com}
#' @seealso \code{\link{lm4d}}, \code{\link{summary.lm.mp}}
#' @examples
#' 
#' # Please see example for lm4d
#' @export
#' 
lm.mp <- function(Y,formula, store.fitted=FALSE) {
    ## Y = eval(formula[[2]], parent.frame())
    X = model.matrix(formula)
    n = dim(X)[1]
    p = dim(X)[2]
    XtX.inv = solve(crossprod(X))
    I.H = diag(n) - X %*% tcrossprod(XtX.inv, X)
    coef = XtX.inv %*% crossprod(X, Y)
    sigma2 = apply(I.H %*% Y, 2, crossprod) / (n-p)
    se.coef = sqrt(diag(XtX.inv) %o% sigma2)
    fitted=if (store.fitted) X %*% coef else NULL
    otpt = list(coef=coef, sigma2=sigma2, se.coef=se.coef, X=X, fitted=fitted) 
    class(otpt) = "lm.mp"
    otpt
}

