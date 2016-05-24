#' @title Cornish-Fisher adjusted variance-covariance VaR
#' 
#' @description Estimates the variance-covariance VaR of a multi-asset portfolio using the Cornish-Fisher adjustment for portfolio-return non-normality, for specified confidence level and holding period.
#' 
#' @param vc.matrix Assumed variance covariance matrix for returns
#' @param mu Vector of expected position returns
#' @param skew Portfolio return skewness
#' @param kurtosis Portfolio return kurtosis
#' @param positions Vector of positions
#' @param cl Confidence level and is scalar or vector
#' @param hp Holding period and is scalar or vector
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' 
#' @examples
#' 
#'    # Variance-covariance for randomly generated portfolio
#'    vc.matrix <- matrix(rnorm(16),4,4)
#'    mu <- rnorm(4)
#'    skew <- .5
#'    kurtosis <- 1.2
#'    positions <- c(5,2,6,10)
#'    cl <- .95
#'    hp <- 280
#'    AdjustedVarianceCovarianceVaR(vc.matrix, mu, skew, kurtosis, positions, cl, hp)
#'    
#' @export
AdjustedVarianceCovarianceVaR <- function(vc.matrix, mu, skew, kurtosis, positions, cl, hp){
  
  # Check that confidence level is read as a row vector
  cl <- as.matrix(cl)
  if (dim(cl)[1] > dim(cl)[2]){
    cl <- t(cl)
  }
  
  # Check that hp is read as a column vector
  hp <- as.matrix(hp)
  if (dim(hp)[1] < dim(hp)[2]){
    hp <- t(hp)
  }
  
  # Check that positions vector read as a scalar or row vector
  positions <- as.matrix(positions)
  if (dim(positions)[1] > dim(positions)[2]){
    positions <- t(positions)
  }
  
  # Check that expected returns vector is read as a scalar or row vector
  mu <- as.matrix(mu)
  if (dim(mu)[1] > dim(mu)[2]){
    mu <- t(mu)
  }
  
  # Check that dimensions are correct
  if (max(dim(mu)) != max(dim(positions))){
    stop("Positions vector and expected returns vector must have same size")
  }
  vc.matrix <- as.matrix(vc.matrix)
  if (max(dim(vc.matrix)) != max(dim(positions))){
    stop("Positions vector and expected returns vector must have same size")
  }
  
  # Check that inputs obey sign and value restrictions
  if (cl >= 1){
    stop("Confidence level must be less than 1")
  }
  if (cl <= 0){
    stop("Confidence level must be greater than 0");
  }
  if (hp <= 0){
    stop("Holding period must be greater than 0");
  }
  
  # Portfolio return standard deviation
  sigma <- positions %*% vc.matrix %*% t(positions)/(sum(positions)^2) # standard deviation of portfolio returns
  # VaR estimation
  z <- double(length(cl))
  adjustment <- z
  VaR <- matrix(0, length(cl), length(hp))
  
  for (i in 1:length(cl)) {
    # Cornish-Fisher adjustment
    z[i] <- qnorm(1 - cl[i], 0, 1)
    adjustment[i] <- (1 / 6) * (z[i] ^ 2 - 1) * skew + (1 / 24) * (z[i] ^ 3 - 3 * z[i]) * (kurtosis - 3) - (1 / 36) * (2 * z[i] ^ 3 - 5 * z[i]) * skew ^ 2
    for (j in 1:length(hp)) {
      VaR[i, j] <- - mu %*% t(positions) * hp[j] - (z[i] + adjustment[i]) * sigma * (sum(positions)^2) * sqrt(hp[j])
    }
  }
  y <- t(VaR)
  return(y)
  
}