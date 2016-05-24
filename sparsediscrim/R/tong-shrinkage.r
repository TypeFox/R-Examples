#' Tong et al. (2012)'s Lindley-type Shrunken Mean Estimator
#'
#' An implementation of the Lindley-type shrunken mean estimator utilized in
#' shrinkage-mean-based diagonal linear discriminant analysis (SmDLDA).
#'
#' @export
#' @references Tong, T., Chen, L., and Zhao, H. (2012), "Improved Mean
#' Estimation and Its Application to Diagonal Discriminant Analysis,"
#' Bioinformatics, 28, 4, 531-537.
#' \url{http://bioinformatics.oxfordjournals.org/content/28/4/531.long}
#' @param x a matrix with \code{n} rows and \code{p} columns.
#' @param r_opt the shrinkage coefficient. If \code{NULL} (default), we calculate
#' the shrinkage coefficient with the formula given just above Equation 5 on page
#' 533 and denoted by \eqn{\hat{r}_{opt}}. We allow the user to specify an
#' alternative value to investigate better approximations.
#' @return vector of length \code{p} with the shrunken mean estimator
tong_mean_shrinkage <- function(x, r_opt = NULL) {
  n <- nrow(x)
  p <- ncol(x)

  # Here, we calculate the approximate "optimal" shrinkage coefficient, r.
  # The formula is given just above Equation 5 and is denoted \hat{r}_{opt}.
  if (is.null(r_opt)) {
    r_opt <- (n - 1) * (p - 2) / n / (n - 3)
  } else {
    r_opt <- as.numeric(r_opt)
  }

  # The sample means of each feature vector.
  xbar <- colMeans(x)
  
  # Tong et al. calculate the mean of the entire matrix, x.
  grand_mean <- mean(x)
  
  # The authors then center the sample mean for each feature vector.
  centered_xbars <- xbar - grand_mean

  # The MLE of the covariance matrix under the assumpton of a multivariate
  # normal population with a diagonal covariance matrix.
  diag_S <- (n - 1) / n * apply(x, 2, var)

  # The term in Equation (6) denoted by:
  # || \bar{x} - \bar{x}_{\dot\dot} ||^2_S
  shrinkage_norm <- sum(centered_xbars^2 / diag_S)

  # Finally, we calculate the shrunken mean given in Equation 6.
  grand_mean + (1 - r_opt / shrinkage_norm) * centered_xbars
}

  
  
