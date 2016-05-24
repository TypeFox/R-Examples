#' @title Hotspots for normal ES
#' 
#' @description Estimates the ES hotspots (or vector of incremental ESs) for a 
#' portfolio assuming individual asset returns are normally distributed, for 
#' specified confidence level and holding period.
#' 
#' @param vc.matrix Variance covariance matrix for returns
#' @param mu Vector of expected position returns
#' @param positions Vector of positions
#' @param cl Confidence level and is scalar
#' @param hp Holding period and is scalar
#' @return Hotspots for normal ES
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' 
#' @examples
#' 
#'    # Hotspots for ES for randomly generated portfolio
#'    vc.matrix <- matrix(rnorm(16),4,4)
#'    mu <- rnorm(4,.08,.04)
#'    positions <- c(5,2,6,10)
#'    cl <- .95
#'    hp <- 280
#'    NormalESHotspots(vc.matrix, mu, positions, cl, hp)
#' 
#' @export
NormalESHotspots <- function(vc.matrix, mu, positions,
                                    cl, hp){
  
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
  if (max(dim(vc.matrix)) != max(dim(positions))){
    stop("Positions vector and variance-covariance matrix must have compatible dimensions")
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
  
  VaR <- - mu %*% t(positions) * hp - qnorm(1 - cl, 0, 1) * 
    (positions %*% vc.matrix %*% t(positions)) * sqrt(hp) # VaR
  n <- 1000 # Number of slives into which tail is divided
  cl0 <- cl # Initial confidence level
  term <- VaR
  delta.cl <- (1 - cl) / n # Increment to confidence level
  for (k in 1:(n - 1)) {
    cl <- cl0 + k * delta.cl # Revised cl
    term <- term - mu %*% t(positions) * hp - qnorm(1 - cl, 0, 1) * 
      (positions %*% vc.matrix %*% t(positions)) * sqrt(hp)
  }
  portfolio.ES <- term/n
  
  # Portfolio ES
  es <- double(length(positions))
  ies <- double(length(positions))
  for (j in 1:length(positions)) {
    x <- positions
    x[j] <- 0
    term[j] <- - mu %*% t(x) * hp - qnorm(1-cl, 0, 1) * x %*% 
      vc.matrix %*% t(x) * sqrt(hp)
    
    for (k in 1:(n - 1)){
      cl <- cl0 + k * delta.cl # Revised cl
      term[j] <- term[j] - mu %*% t(x) * hp - qnorm(1-cl, 0, 1) * x %*% 
        vc.matrix %*% t(x) * sqrt(hp)
    }
    es[j] <- term[j]/n # ES on portfolio minus position j
    ies [j] <- portfolio.ES - es[j] # Incremental ES
    
  }
  y <- ies
  return(ies)
}
