#' F-tests for massively parallel linear models
#' 
#' Performs F-tests for removing one or more terms from each of a large number
#' of models with common design matrix.
#' 
#' 
#' @param formula a formula such as "\code{Y ~ X}", where \code{Y} is an \eqn{n
#' \times V} response matrix and \code{X} is an \eqn{n \times p} design matrix
#' common to all \eqn{V} models.
#' @param which number or vector indicating which column(s) of the model matrix
#' are to be tested for removal from the model.
#' @return \item{F}{F-statistics for each of the models.} \item{df1}{numerator
#' degrees of freedom.} \item{df2}{denominator degrees of freedom.}
#' \item{pvalue}{upper-tailed p-value.} \item{X}{design matrix.}
#' @author Philip Reiss \email{phil.reiss@@nyumc.org} and Lei Huang
#' \email{huangracer@@gmail.com}
#' @seealso \code{\link{lm.mp}}, \code{\link{permF.mp}}
#' @examples
#' 
#' Y = matrix(rnorm(6000), nrow=20)
#' X = rnorm(20)
#' t2 = F.mp(Y~X, which=2)
#' @export
F.mp <-
function(formula, which)    {
    Y = eval(formula[[2]], parent.frame())
    X = model.matrix(formula)
    n = dim(X)[1]
    p = dim(X)[2]
    df1 = length(which)
    XtX.inv = solve(crossprod(X))
    I.H = diag(n) - X %*% tcrossprod(XtX.inv, X)
    rss1 = apply(I.H %*% Y, 2, crossprod)
    X0 = matrix(X[ , -which], nrow = nrow(X))
    XtX.inv0 = solve(crossprod(X0))
    I.H0 = diag(n) - X0 %*% tcrossprod(XtX.inv0, X0)
    rss0 = apply(I.H0 %*% Y, 2, crossprod)
    F = ((n-p)/df1) * (rss0-rss1) / rss1
    otpt = list(F=F, df1=df1, df2=n-p, pvalue=pf(F, df1, n-p, lower.tail=FALSE), X=X)
    otpt
}

