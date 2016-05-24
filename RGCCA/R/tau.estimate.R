#' Estimation of the optimal shrinkage parameters as described in [1,2] and implemented in \code{\link[SHIP:SHIP-package]{SHIP}} [2]. 
#' @param x  Data set on which the covariance matrix is estimated.
#' @return \item{tau}{Optimal shrinkage intensity parameter}
#' @title Optimal shrinkage intensity parameters.
#' @references [1] Schaefer J. and Strimmer K., 2005. A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @references [2] Jelizarow M., Guillemot V., Tenenhaus A., Strimmer K., Boulesteix A.-L., 2010. Over-optimism in bioinformatics: an illustration. Bioinformatics 26:1990-1998.
#' @export tau.estimate

tau.estimate <- function(x) {
  if (is.matrix(x) == TRUE && is.numeric(x) == FALSE) 
    stop("The data matrix must be numeric!")
  p <- NCOL(x)
  n <- NROW(x)
  covm <- cov(x)
  corm <- cor(x)
  xs <- scale(x, center = TRUE, scale = TRUE)
  v <- (n/((n - 1)^3)) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
  diag(v) <- 0
  m <- matrix(rep(apply(xs^2, 2, mean), p), p, p)
  I <- diag(NCOL(x))
  d <- (corm - I)^2
  tau <- (sum(v))/sum(d)
  tau <- max(min(tau, 1), 0)
  return(tau)
 }

