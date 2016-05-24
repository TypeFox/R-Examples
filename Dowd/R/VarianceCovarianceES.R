#' @title Variance-covariance ES for normally distributed returns
#' 
#' @description Estimates the variance-covariance VaR of a
#' portfolio assuming individual asset returns are normally distributed, 
#' for specified confidence level and holding period.
#' 
#' @param vc.matrix Variance covariance matrix for returns
#' @param mu Vector of expected position returns
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
#'    # Variance-covariance ES for randomly generated portfolio
#'    vc.matrix <- matrix(rnorm(16), 4, 4)
#'    mu <- rnorm(4)
#'    positions <- c(5, 2, 6, 10)
#'    cl <- .95
#'    hp <- 280
#'    VarianceCovarianceES(vc.matrix, mu, positions, cl, hp)
#' 
#' @export
VarianceCovarianceES <- function(vc.matrix, mu, positions, cl, hp){
  
  # Check that cl is read as a row vector
  cl <- as.matrix(cl)
  if (dim(cl)[1] > dim(cl)[2]) {
    cl <- t(cl)
  }
  
  # Check that hp is read as a column vector
  hp <- as.matrix(hp)
  if (dim(hp)[1] < dim(hp)[2]) {
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
  
  # VaR and ES estimation
  VaR <- matrix(0, length(cl), length(hp))
  term <- matrix(0, length(cl), length(hp))
  es <- matrix(0, length(cl), length(hp))
  cl0 <- double(length(cl))
  delta.cl <- double(length(cl))
  for (i in 1:length(cl)) {
    for (j in 1:length(hp)) {
      VaR[i,j] <- - mu %*% t(positions) * hp[j] - qnorm(1-cl[i],0,1) * 
        (positions%*%vc.matrix%*%t(positions))*sqrt(hp[j]) # VaR
      # ES Estimation
      n <- 1000 # Number of slives into which tail is divided
      cl0[i] <- cl[i] # Initial confidence level
      term[i, j] <- VaR[i, j]
      delta.cl[i] <- (1 - cl[i]) / n # Increment to confidence level as each 
      # slice is taken
      
      for (k in 1:(n - 1)) {
        
        cl[i] <- cl0[i] + k * delta.cl[i] # Revised cl
        term[i, j] <- term[i, j] - mu %*% t(positions) * hp[j] - 
          (qnorm(1-cl[i],0,1)) * (positions%*%vc.matrix%*%t(positions))*sqrt(hp[j])
      }
      es[i, j] <- term[i, j]/n
      
    }
  }
  y <- t(es)
  return(y)
  
}
