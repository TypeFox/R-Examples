predict.multi_tpfit <-
function(object, lags, byrow=TRUE, ...) {
  # Transition probabilities prediction for spatial Markov Chain in nD
  #
  #    lags vector or matrix of lags
  #  object list with a matrix of estimated transition rates

  res <- list()
  class(res) <- "multi_transiogram"
  res$type <- "Theoretical"

  nc <- length(object$coefficients)
  if (is.vector(lags)) {
    if (nc != length(lags)) stop ("wrong the length of \"lags\"")
    lags <- matrix(lags, nrow = nc)
  }
  else {
    if (!is.matrix(lags)) stop ("\"lags\" must be a matrix or a vector")
    lags <- if (byrow) { t(lags) } else { lags }
    if (dim(lags)[1] != nc) stop ("wrong number of ", ifelse(byrow, "columns", "rows"))
  }
  storage.mode(lags) <- "double"

  nk <- length(object$prop)
  nr <- dim(lags)[2]

  coefficients <- unlist(object$coefficients)
  res$Tmat <- array(0, dim = c(nk, nk, nr))

  if (!is.null(object$rotation)) {
    dire.mat <- diag(, nc)
    dire.mat <- .C('rotaxes', nc = as.integer(nc), ang = as.double(object$rotation), 
                   res = as.double(dire.mat), PACKAGE = "spMC")$res
    dire.mat <- t(matrix(dire.mat, nc, nc))
    rotlags <- .C('fastMatProd', nr = as.integer(nc), ni = as.integer(nc),
                  mat1 = as.double(dire.mat), nc = as.integer(nr), mat2 = as.double(lags),
                  res = as.double(vector("numeric", prod(dim(lags)))),
                  PACKAGE = "spMC")$res
    res$Tmat <- .C('predMULTI', coefficients = as.double(coefficients),
                   prop = as.double(object$prop), lags = as.double(rotlags),
                   nk = as.integer(nk), nc = as.integer(nc), nr = as.integer(nr),
                   mypred = as.double(res$Tmat), PACKAGE = "spMC")$mypred
  }
  else {
    res$Tmat <- .C('predMULTI', coefficients = as.double(coefficients),
                   prop = as.double(object$prop), lags = as.double(lags),
                   nk = as.integer(nk), nc = as.integer(nc), nr = as.integer(nr),
                   mypred = as.double(res$Tmat), PACKAGE = "spMC")$mypred
  }

  res$Tmat <- array(res$Tmat, dim = c(nk, nk, nr))
  colnames(res$Tmat) <- rownames(res$Tmat) <- names(object$prop)
  res$lags <- t(lags)
  return(res)
}
