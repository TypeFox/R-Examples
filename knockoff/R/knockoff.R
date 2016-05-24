#' The knockoff filter for controlling the false discovery rate
#'
#' @docType package
#' @name knockoff
#' @import RJSONIO
NULL

#' Knockoff filter
#' 
#' This function runs the knockoff procedure from start to finish, creating 
#' the knockoffs, computing the test statistics, and selecting variables. 
#' It is the main entry point for the knockoff package.
#' 
#' @param X matrix or data frame of predictors
#' @param y response vector
#' @param fdr target FDR (false discovery rate)
#' @param statistic the test statistic (by default, a lasso statistic). See the
#'  Details section for more information.
#' @param threshold either 'knockoff' or 'knockoff+'.
#' @param knockoffs either equicorrelated knockoffs ('equicorrelated') or 
#'  knockoffs optimized using semidefinite programming ('sdp')
#' @param normalize whether to scale the data columns to have unit norm. Only
#'  disable this if your data is already normalized.
#' @param randomize whether randomization is to be used when constructing
#'  knockoffs and (when p < n < 2p) augmenting the model with extra rows
#'
#' @return An object of class "knockoff.result". This object is a list 
#'  containing at least the following components:
#'  \item{knockoff}{matrix of knockoff variables}
#'  \item{statistic}{computed test statistic}
#'  \item{threshold}{computed selection threshold}
#'  \item{selected}{named vector of selected variables}
#'
#' @details The default test statistic is \code{knockoff.stat.lasso_signed_max}.
#' Other useful test statistics include \code{knockoff.stat.fs} and 
#' \code{knockoff.stat.fs_omp}. It is also possible to provide your own test 
#' statistic (for an example, see the vignette).
#' 
#' To use SDP knockoffs, you must have a Python installation with 
#' CVXPY. For more information, see the vignette on SDP knockoffs:
#' 
#' \code{vignette('sdp', package='knockoff')}
#' 
#' @export
knockoff.filter <- function(X, y, fdr=0.20, statistic=NULL, 
                            threshold=c('knockoff','knockoff+'),
                            knockoffs=c('equicorrelated','sdp'),
                            normalize=TRUE, randomize=FALSE) {
  # Parameter defaults.
  if (is.null(statistic)) statistic = knockoff.stat.lasso_signed_max
  
  # Validate input types.
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)  
  } else if (is.matrix(X))
    X.names = colnames(X)
  else
    stop('Input X must be a matrix or data frame')
  y = as.vector(y)
  
  # Validate input dimensions.
  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n)
  if (n <= p)
    stop('Input X must have dimensions n > p')
  else if (n < 2*p) {
    warning('Input X has dimensions p < n < 2p. ',
            'Augmenting the model with extra rows.')
    X.svd = svd(X, nu=n, nv=0)
    u2 = X.svd$u[,(p+1):n]
    sigma = sqrt(mean((t(u2) %*% y)^2)) # = sqrt(RSS/(n-p))
    X = rbind(X, matrix(0, 2*p-n, p))
    if (randomize)
      y.extra = rnorm(2*p-n, sd=sigma)
    else
      y.extra = with_seed(0, rnorm(2*p-n, sd=sigma))
    y = c(y, y.extra)
  }
  
  # Normalize X, if necessary.
  if (normalize)
    X = normc(X)
  
  # Run the knockoff filter!
  X.knockoff = knockoff.create(X, method=knockoffs, randomize=randomize)
  W = statistic(X, X.knockoff, y)
  t = knockoff.threshold(W, fdr, threshold)
  
  selected = which(W >= t)
  if (!is.null(X.names))
    names(selected) = X.names[selected]
  
  # Package up the results.
  structure(list(call = match.call(),
                 knockoff = X.knockoff,
                 statistic = W,
                 threshold = t,
                 selected = selected),
            class = 'knockoff.result')
}

#' Threshold for the knockoff filter
#' 
#' Computes the threshold for the knockoff filter.
#' 
#' @param W the test statistics
#' @param fdr target FDR (false discovery rate)
#' @param method either 'knockoff or 'knockoff+'
#' @return The threshold for variable selection.
#' 
#' @export
knockoff.threshold <- function(W, fdr=0.20, method=c('knockoff','knockoff+')) {
  offset = switch(match.arg(method),
                  knockoff = 0, `knockoff+` = 1)
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

#' @export
#' @keywords internal
print.knockoff.result <- function(x, ...) {
  cat('Call:\n')
  print(x$call)
  cat('\nSelected variables:\n')
  print(x$selected)
}