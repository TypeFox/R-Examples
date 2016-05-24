#' For given series of operators \eqn{(A_k)_{k \in D}} where \eqn{D = }\code{A$lags}
#' Computes operator norms for each element and returns a series \eqn{\|A_k\|}
#'
#' @title Compute operator norms of elements of a filter
#' @param A the series of operators
#' @param type matrix norm to be used as in \link[base]{norm}
#' @return List with \code{lags} and \code{norms} vetors
#' @export
timedom.norms = function(A, type='2'){
  if (!is.timedom(A))
    stop ("A must be a time domain filter")
  
  mynorm = norm
  if (type=='2')
    mynorm = norm.spec
  R = list()
  R$lags = A$lags
  v = c()
  D = dim(A$operators)
  for (i in 1:D[1])
    v = c(v,mynorm(matrix(A$operators[i,,],D[2],D[3])))
  R$norms = v
  R
}

