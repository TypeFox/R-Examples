#' Estimates VaR with Box-Cox transformation
#' 
#' Function estimates the VaR of a portfolio assuming P and L data set transformed
#' using the BoxCox transformation to make it as near normal as possible, for 
#' specified confidence level and holding period implied by data frequency.
#' 
#' @param PandLdata Daily Profit/Loss data
#' @param cl Confidence Level. It can be a scalar or a vector.
#' @return Estimated Box-Cox VaR. Its dimension is same as that of cl
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' Hamilton, S. A. and Taylor, M. G. A Comparision of the Box-Cox 
#' transformation method and nonparametric methods for estimating quantiles
#' in clinical data with repeated measures. J. Statist. Comput. Simul., vol. 
#' 45, 1993, pp. 185 - 201.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates Box-Cox VaR
#'    a<-rnorm(100)
#'    BoxCoxVaR(a,.95)
#' 
#' @import forecast
#' 
#' @export
BoxCoxVaR <- function(PandLdata, cl){
  # Check that inputs have correct dimensions
  cl <- as.matrix(cl)
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  if (min(cl.row, cl.col) > 1) {
    stop("Confidence level must be a scalar or a vector")
  }
  
  if (cl.row > cl.col) {
    cl <- t(cl)
  }
  
  # Check that inputs obey sign and value restrictions
  if (max(cl) >= 1){
    stop("Confidence level(s) must be less than 1")
  }
  if (min(cl) <= 0){
    stop("Confidence level(s) must be greater than 0")
  }
  # Transform data and obtain lambda
  loss.data <- -PandLdata
  loss.data <- loss.data - min(loss.data) + 1
  loss.data <- sort(loss.data)
  lambda <- BoxCox.lambda(loss.data, method="loglik")
  transdat <- BoxCox(loss.data, lambda)
  
  # Alternative method: 
  # for dependence only on MASS and not on forecast package (not working yet!)
  # model <- lm(loss.data~1)
  # boxcox <- boxcox(model,plotit=FALSE)
  # lambda <- with(bc, x[which.max(y)])
  # box cox transformation
  # if(lambda == 0) {
  #  transdat <- log(loss.data)
  # } else {
  #  transdat <- (loss.data^lambda -1)/lambda
  # }
  
  # Estimate mean and standard deviation of transformed data
  mu <- mean(transdat)
  sigma <- sd(transdat)
  VaR <- double(length(cl))
  for(i in 1:length(cl)){
    VaR[i] <- (1 + lambda * (mu + sigma * qnorm(cl[i]))) ^ (1 / lambda) + min(-PandLdata) - 1 # i-th VaR
  }
  return(VaR)
}