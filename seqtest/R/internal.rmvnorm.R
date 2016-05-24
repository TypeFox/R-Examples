##########################################################################################################
#
# seqtest: Sequential Triangular Test
#
# Internal function: rmvnorm from the package "mvtnorm"
#
# Author: Friedrich Leisch
#
internal.rmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                             method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE) {

  #-----------------------------------------------------------------------------------
  # Input Check

  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {

    stop("sigma must be a symmetric matrix")

  }

  ###

  if (length(mean) != nrow(sigma)) {

    stop("mean and sigma have non-conforming size")

  }

  #-----------------------------------------------------------------------------------
  # Main function

  method <- match.arg(method)

  if (method == "eigen") {

    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {

      warning("sigma is numerically not positive definite")

    }

    R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))

  }

  ###

  if (method == "svd") {

    s. <- svd(sigma)

    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {

      warning("sigma is numerically not positive definite")

    }

    R <- t(s.$v %*% (t(s.$u) * sqrt(s.$d)))

  }

  ###

  if (method == "chol") {

    R <- chol(sigma, pivot = TRUE)
    R <- R[, order(attr(R, "pivot"))]

  }

  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
  retval <- sweep(retval, 2, mean, "+")

  colnames(retval) <- names(mean)

  return(retval)

}
