#---------------------------------------------------------------------
# l2boost data simulations
#
# This is to encapsulate all simulation dataset in the application file
#---------------------------------------------------------------------

#' @title multivariate normal data simulations.
#'
#' @description Create simulated dataset from a multivariate normal. Used to recreate data simulations
#'  from Ehrlinger and Ishwaran (2012).
#' 
#' @param n number of observations
#' @param p number of coordinate directions in the design matrix
#' @param beta a "true" beta vector of length p (default=NULL) See details.
#' @param which.beta indicator vector for which beta coefficients to include as signal in simulation (default=NULL) see details
#' @param rho correlation coefficient between coordinate directions 
#'
#' @details By default, mvnorm.l2boost creates a data set of n multivariate normal random observations of p covariates 
#' (see MASS:mvrnorm). The correlation matrix is constructed with 1 on the diagonals and the correlation 
#' coefficient \emph{rho} on the off diagonals. 
#' 
#' The response is constructed as follows: If a true beta vector is not supplied, the first 10 beta coefficients carry 
#' the signal with a value of 5, and the remaining p-10 values are set to zero. Given a \emph{beta.true} vector, all 
#' values are used as specified. The coefficent vector is truncated to have \emph{p} signal terms if 
#' length(\emph{beta.true}) > \emph{p}, and noise coordinates are added if length(\emph{beta.true}) < \emph{p}.
#' 
#' It is possible to pass an indicator vector \emph{which.beta} to select specific signal elements from the 
#' full vector \emph{beta.true}. 
#' 
#' @return 
#' \itemize{
#' \item call Matched function call
#' \item x design matrix of size \emph{n} x \emph{p}
#' \item y response vector of length \emph{n}
#' }
#' @examples
#' #--------------------------------------------------------------------------
#' # Example: Multivariate normal data simulation
#'
#' # Create a (reproducable) data set of size 100 x 100
#' set.seed(1024)
#' n<- 100
#' p<- 100
#' 
#' # Set 10 signal variables using a uniform beta=5, the remaining (p-10)=90 are
#' # set to zero indicating random noise.  
#' beta <- c(rep(5,10), rep(0,p-10))
#' 
#' # Example with orthogonal design matrix columns (orthogonal + noise)
#' ortho.data <- mvnorm.l2boost(n, p, beta)
#' cbind(ortho.data$y[1:5],ortho.data$x[1:5,])
#' 
#' # Example with correlation between design matrix columns
#' corr.data <- mvnorm.l2boost(n, p, beta, rho=0.65)
#' cbind(corr.data$y[1:5],corr.data$x[1:5,])
#' 
#' 
#' @references Ehrlinger J., and Ishwaran H. (2012). "Characterizing l2boosting" \emph{Ann. Statist.}, 40 (2), 1074-1101
#' 
#' @export mvnorm.l2boost
#' @importFrom MASS mvrnorm
mvnorm.l2boost <- function(n=100,p=100, beta=NULL, which.beta=NULL, rho=0){
  call<-match.call()
  if(rho == 0){
    sigma <- diag(1, p)
  }else{
    sigma <- matrix(rho, p, p)
    diag(sigma) <- 1
  }
  x <- mvrnorm(n, mu = rep(0, p), Sigma = sigma)
  
  beta.true <- rep(0, p)   
  if(is.null(beta)){
    if(is.null(which.beta)) beta.true[1:10] <- 5 else beta.true[which.beta] <- 5 
  }else if(length(beta) == p){
    beta.true = beta
  }else if(is.null(which.beta)){
    beta.true[1:length(beta)] <- beta
  }else{
    beta.true[which.beta] <- beta
  }
  y <- x %*% beta.true + rnorm(n)
  return(list(x=x, y=y, call =call))
}
