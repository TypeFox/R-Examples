#' @title Hotspots for VaR adjusted by Cornish-Fisher correction
#' 
#' @description Estimates the VaR hotspots (or vector of incremental VaRs) for a
#' portfolio with portfolio return adjusted for non-normality by Cornish-Fisher 
#' corerction, for specified confidence level and holding period.
#' 
#' @param vc.matrix Variance covariance matrix for returns
#' @param mu Vector of expected position returns
#' @param skew Return skew
#' @param kurtosis Return kurtosis
#' @param positions Vector of positions
#' @param cl Confidence level and is scalar
#' @param hp Holding period and is scalar
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' 
#' @examples
#' 
#'    # Hotspots for ES for randomly generated portfolio
#'    vc.matrix <- matrix(rnorm(16),4,4)
#'    mu <- rnorm(4)
#'    skew <- .5
#'    kurtosis <- 1.2
#'    positions <- c(5,2,6,10)
#'    cl <- .95
#'    hp <- 280
#'    AdjustedNormalVaRHotspots(vc.matrix, mu, skew, kurtosis, positions, cl, hp)
#'    
#' @export
AdjustedNormalVaRHotspots <- function(vc.matrix, mu, skew, kurtosis, positions, cl, hp){
  
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
    stop("Positions vector and variance-covariance matrix must have compatible dimensions.")
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
  
  # VaR and ES estimation
  z <- qnorm(1 - cl, 0 ,1)
  sigma <- positions %*% vc.matrix %*% t(positions)/(sum(positions)^2) # Initial standard deviation of portfolio returns
  adjustment <- (1 / 6) * (z ^ 2 - 1) * skew + (1 / 24) * (z ^ 3 - 3 * z) * (kurtosis - 3) - (1 / 36) * (2 * z ^ 3 - 5 * z) * skew ^ 2
  VaR <- - mu %*% t(positions) * hp - (z + adjustment) * sigma * (sum(positions)^2) * sqrt(hp)
  
  # VaR
  x <- double(length(positions))
  sigma <- double(length(positions))
  iVaR <- double(length(positions))
  for (i in 1:length(positions)){
    x <- positions
    x[i] <- 0
    sigma[i] <- x %*% vc.matrix %*% t(x)/sum(x)^2 # standard deviation of portfolio returns
    iVaR[i] <- VaR + mu %*% t(x) %*% hp + (z + adjustment) * sigma[i] * (sum(x))^2 * sqrt(hp) # Incremental VaR
  }
  y <- iVaR
  return(y)
  
}