#' Estimates ES with Box-Cox transformation
#' 
#' Function estimates the ES of a portfolio assuming P and L data set transformed
#' using the BoxCox transformation to make it as near normal as possible, for 
#' specified confidence level and holding period implied by data frequency.
#' 
#' @param loss.data Daily Profit/Loss data
#' @param cl Confidence Level. It can be a scalar or a vector.
#' @return Estimated Box-Cox ES. Its dimension is same as that of cl
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
#'    a<-rnorm(200)
#'    BoxCoxES(a,.95)
#' 
#' @import forecast
#' 
#' @export
BoxCoxES <- function(loss.data, cl){
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
  # ES Estimation
  VaR <- BoxCoxVaR(loss.data, cl) # HS VaR
  k <- which(VaR<loss.data) # Finds indices of tail loss data
  tail.losses <- loss.data[k] # Creates data set of tail loss observations
  ES <- mean(tail.losses) # ES
  return(ES) 
}