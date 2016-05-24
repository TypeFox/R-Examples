# Compute crude consistent estimate ,that's, equation (6) in the reference paper.
# the initial sample size m1 is adjusted from 40 to 60.
y_crude_estimator <- function(mk_correlated_standard_norm,
                              invcdfnames, paramslists) {
  if (class(mk_correlated_standard_norm) != "matrix") {
    stop("mk_correlated_standard_norm must be matrix !")
  }
  if (length(invcdfnames) != length(paramslists)) {
    stop("inversecdfs should have the same length paramslists !")
  }
  ndim <- ncol(mk_correlated_standard_norm)
  mk <- nrow(mk_correlated_standard_norm)
  if (mk >= 60) {
    stop("the number of observations must be less than 60!")
  }
  transform_mat <- NULL
  for (i in 1:ndim) {

  funcall <- as.call(c(as.name(invcdfnames[i]),
                      list(pnorm(mk_correlated_standard_norm)[ ,i]), paramslists[[i]]))
  transform_mat <- cbind(transform_mat, eval(funcall))
  }

  res <- matrix(rep(0, ndim * ndim), nrow = ndim)
  diag(res) <- 1
  for (i in 1:(ndim-1))
    for (j in (i+1):ndim)
      if (length(which(!duplicated(transform_mat[ ,i])[-1]))==0 ||
           length(which(!duplicated(transform_mat[ ,j])[-1]))==0)
      res[j,i] <- res[i,j] <- cor(transform_mat[ ,i], transform_mat[ ,j])
      else
      res[j,i] <- res[i,j] <- 0
  res
}
