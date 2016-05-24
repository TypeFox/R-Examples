#' Performs nfolds-cross validation
#'
#' This function can be used to select a value of lam that performs well according to a user-specified measure of error.
#' 
#' @param pathObj output of hierband.path
#' @param x \code{n} by \code{p} matrix
#' @param errfun a user-specified function measuring the loss incurred by estimating \code{est} (first argument)
#' when the true covariance matrix is \code{true} (second argument).  Default: Squared Frobenius norm.
#' @param nfolds number of folds (default: 5)
#'
#' @return 
#' \describe{
#' \item{errs: }{A \code{nlam}-by-\code{nfolds} matrix of errors.  \code{errs[i,j]} is error incurred in using \code{lamlist[i]} on fold \code{j}}
#' \item{m: }{CV error error for each value of lambda.}
#' \item{se: }{Standard error (estimated over folds) for each value of lambda}
#' \item{lam.best: }{Value of \code{lamlist} minimizing CV error.}
#' \item{ibest: }{Index of \code{lamlist} minimizing CV error.}
#' \item{lam.1se.rule: }{Selected value of lambda using the one-standard-error rule, a common heuristic that favors a sparser model when there isn't strong evidence against it.}
#' \item{i.1se.rule: }{Index of \code{lamlist} of one-standard-error rule.}
#' }
#'
#' @export
#'
#' @seealso \code{\link{hierband}} \code{\link{hierband.path}}
#' @examples
#' set.seed(123)
#' p <- 100
#' n <- 50
#' K <- 10
#' true <- ma(p, K)
#' x <- matrix(rnorm(n*p), n, p) %*% true$A
#' Sighat <- cov(x)
#' path <- hierband.path(Sighat)
#' cv <- hierband.cv(path, x)
#' fit <- hierband(Sighat, lam=cv$lam.best)
#' \dontrun{
#' plot(cv$m, main="CV Frob Error", type="b")
#' lines(cv$m+cv$se, main="CV Frob Error")
#' lines(cv$m-cv$se, main="CV Frob Error")
#' abline(v=c(cv$ibest,cv$i.1se.rule), lty=2)
#' }
hierband.cv <- function(pathObj, x, errfun=function(est,true) sum((est-true)^2), nfolds=5) {
  n <- nrow(x)
  folds <- MakeFolds(n, nfolds)
  nlam <- length(pathObj$lamlist)
  errs <- matrix(NA, nlam, nfolds)
  for (i in seq(nfolds)) {
    # train on all but i-th fold (and use settings from pathObj):
    fit <- hierband.path(cov(x[-folds[[i]],]), lamlist=pathObj$lamlist, w=pathObj$w, delta=pathObj$delta)
    # evaluate this on left-out fold:
    Sig.te <- cov(x[folds[[i]],])
    for (j in seq(nlam)) errs[j, i] <- errfun(fit$P[,,j], Sig.te)
  }
  m <- rowMeans(errs)
  se <- apply(errs,1,sd)/sqrt(nfolds)
  ibest <- which.min(m)
  i.1se.rule <- min(which(m < m[ibest]+se[ibest]))
  list(errs=errs, m=m, se=se, lam.best=fit$lamlist[ibest], ibest=ibest,
       lam.1se.rule=fit$lamlist[i.1se.rule], i.1se.rule=i.1se.rule)
}

#' Make folds for cross validation
#' 
#' Divides the indices \code{1:n} into \code{nfolds} random folds of about the same size.
#' 
#' @param n sample size
#' @param nfolds number of folds
MakeFolds <- function(n, nfolds) {
  nn <- round(n / nfolds)
  sizes <- rep(nn, nfolds)
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
  b <- c(0, cumsum(sizes))
  ii <- sample(n)
  folds <- list()
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds
}
