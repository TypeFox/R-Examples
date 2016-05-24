#' Estimates VaR by principal components analysis
#' 
#' Estimates the VaR of a multi position portfolio by principal components analysis, using chosen number of principal components and a specified confidence level or range of confidence levels.
#' 
#' @param Ra Matrix return data set where each row is interpreted as a set of daily observations, and each column as the returns to each position in a portfolio
#' @param position.data Position-size vector, giving amount invested in each position
#' @param number.of.principal.components Chosen number of principal components
#' @param cl Chosen confidence level
#' @return VaR
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes PCA VaR
#'    Ra <- matrix(rnorm(4*6),4,6)
#'    position.data <- rnorm(6)
#'    PCAVaR(Ra, position.data, 2, .95)
#'
#' @export
PCAVaR <- function(Ra, position.data, number.of.principal.components, cl){
  # Check that inputs have correct dimensions
  return.data<-as.matrix(Ra)
  m <- dim(return.data)[1]
  n <- dim(return.data)[2]
  if (min(m, n) == 1) {
    stop("Input data set has insufficient dimensionality")
  }
  if (number.of.principal.components < 0) {
    stop("Number of principal components must be positive")
  }
  if (number.of.principal.components > n) {
    stop("Number of principal components cannot exceed number of positions")
  }
  # Check that dimensions of position data and return data match
  if (n != length(position.data)) {
    stop("Dimensions of return data and position data should be compatible")
  }
  
  # Principal components estimation
  a <- svd(return.data)  # SVD; provides U and V
  S.diag <- sort(a$d, decreasing = TRUE)[1:number.of.principal.components] # Creates diagonal for S matrix
  # Following condition for the fact that scalar argument to diag returns 
  # identity matrix
  if (length(S.diag) == 1) {
    S.diag <- as.matrix(S.diag)
  } else {
    S.diag <- diag(S.diag)
  }
  S.diag <- diag(S.diag)
  S <- matrix(0, min(m, n), min(m, n))
  S[1:number.of.principal.components,1:number.of.principal.components] <- S.diag # Creates S matrix with diagonal S.diag
  synthetic.PandL.data <- a$u %*% S %*% t(a$v) %*% position.data
  
  y <- HSVaR(synthetic.PandL.data, cl)
  return(y)
}