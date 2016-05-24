#' @title Variance-covariance VaR for normally distributed returns
#' 
#' @description Estimates the variance-covariance VaR of a
#' portfolio assuming individual asset returns are normally distributed, 
#' for specified confidence level and holding period.
#' 
#' @param vc.matrix Assumed variance covariance matrix for returns
#' @param mu Vector of expected position returns
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
#'    # Variance-covariance VaR for randomly generated portfolio
#'    vc.matrix <- matrix(rnorm(16),4,4)
#'    mu <- rnorm(4)
#'    positions <- c(5,2,6,10)
#'    cl <- .95
#'    hp <- 280
#'    VarianceCovarianceVaR(vc.matrix, mu, positions, cl, hp)
#'    
#' @seealso AdjustedVarianceCovarianceVaR
#' @export
VarianceCovarianceVaR <- function(vc.matrix, mu, positions, cl, hp){
  
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
  
  # Check that positions vector is read as a scalar or row vector
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
  
  # VaR estimation
  VaR <- matrix(0, length(cl), length(hp))
  for (i in 1:length(cl)) {
    for (j in 1:length(hp)) {
      VaR[i, j] <- - mu %*% t(positions) * hp[j] - qnorm(1-cl[i], 0, 1) * (positions %*% vc.matrix %*% t(positions)) * sqrt(hp[j])
    }
  }
  y <- t(VaR)
  return(y)
}