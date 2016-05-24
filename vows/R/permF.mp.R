#' Permutation F-tests for massively parallel linear models
#' 
#' Performs permutation F-tests for parallel linear models with a common design
#' matrix.  Currently restricted to testing with the intercept-only model as
#' the null hypothesis.  The permutation method controls the familywise error
#' rate (FWER) at a desired level; see Details.
#' 
#' The observed F-statistics are referred to a permutation distribution of the
#' maximum F-statistic over all \eqn{V} tests.  This is a standard approach to FWER
#' control in neuroimaging (Nichols and Holmes, 2001).
#' 
#' @param formula a formula such as "\code{Y ~ X}", where \code{Y} is an \eqn{n
#' \times V} response matrix and \code{X} is an \eqn{n \times p} design matrix
#' common to all \eqn{V} models.
#' @param nperm number of permutations.
#' @param alpha level at which to control the FWER.
#' @param report.every parameter controlling how often to report the number of
#' permutations performed; by default, every \code{50}.
#' @return \item{maxF.perm}{maximal F-statistics obtained from each of the
#' permuted data sets.} \item{F.obs}{the observed F-statistics.}
#' \item{threshold}{critical value obtained from the permutations.}
#' \item{pvalue}{adjusted (familywise error rate-controlling) p-values.}
#' @author Philip Reiss \email{phil.reiss@@nyumc.org} and Lei Huang
#' \email{huangracer@@gmail.com}
#' @seealso \code{\link{F.mp}}
#' @references Nichols, T. E., and Holmes, A. P. (2001).  Nonparametric
#' permutation tests for functional neuroimaging: a primer with examples.
#' \emph{Human Brain Mapping}, 15, 1--25.
#' @examples
#' 
#' Y = matrix(rnorm(6000), nrow=20)
#' X = rnorm(20)
#' t3 = permF.mp(Y~X)
#' @export
permF.mp <-
function(formula, nperm=499, alpha=.05, report.every=50)  {
    Y = eval(formula[[2]], parent.frame())
    X = model.matrix(formula)
    n = dim(X)[1]; p = dim(X)[2]
    I.H = diag(n) - X %*% solve(crossprod(X), t(X))
    rss0 = apply(scale(Y, TRUE, FALSE), 2, crossprod)
    maxF.perm = c()
    for (i in 1:nperm)  {
        if (i %% report.every==0) cat("Permutation", i, "\n")
        Y.perm = Y[sample(1:n), ]
        rss1.perm = apply(I.H %*% Y.perm, 2, crossprod)
        F.perm = ((n-p)/(p-1)) * (rss0-rss1.perm) / rss1.perm
        maxF.perm[i] = max(F.perm)
    }
    rss1 = apply(I.H %*% Y, 2, crossprod)
    F.obs = ((n-p)/(p-1)) * (rss0-rss1) / rss1
    threshold = quantile(maxF.perm, 1-alpha)
    pvalue = (rowSums(outer(F.obs, maxF.perm, function(x,y) x<y)) + 1) / (nperm + 1)
    return(list(maxF.perm=maxF.perm, F.obs=F.obs, threshold=threshold, pvalue=pvalue))
}

