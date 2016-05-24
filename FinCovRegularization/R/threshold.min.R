#' @title Minimum threshold constant
#'
#' @description
#' This function is for determining the minimum constant in the threshold
#' that guarantees the positive definiteness of the estimator.
#'
#' @param sigma a covariance matrix
#' @param method a character, indicating thresholding method: "soft" or "hard"
#' @return minimum constant for thresholding
#' @references "High-Dimensional Covariance Estimation" by Mohsen Pourahmadi
#' @keywords internal
#' @importFrom stats uniroot
#' @export

threshold.min <- function(sigma, method = "hard") {
  if ((method %in% c("hard","soft")) == FALSE) {
    stop("This function only support two thresholding methods: hard and soft")
  }

  mineigen <- function(sigma, threshold, method) {
    if (method == "hard") {
      COV <- hard.thresholding(sigma, threshold)
    } else if (method == "soft") {
      COV <- soft.thresholding(sigma, threshold)
    }
    mineigen <- min(eigen(COV)$values)
    return(mineigen)
  }
  f <- function(x) mineigen(sigma, threshold = x, method)

  if (f(0) * f(max(sigma)) < 0) {
    r <- uniroot(f, c(0, max(sigma)), tol = sqrt(.Machine$double.eps))
    threshold.min <- max(0, r$root)
    return(threshold.min)
  } else {
    threshold.min <- 0
    return(threshold.min)
  }
}
