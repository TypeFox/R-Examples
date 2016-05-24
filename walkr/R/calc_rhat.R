#' R hat 
#' 
#' This function calculates the rhat
#' of each parameter given the list of chains 
#' that walkr produces. Since this is just an 
#' internal function, I'll document it more later.
#' 
#' @param x is the list of chains 
#' @importFrom stats var
#' 
#' @return a vector of rhats 
#' 
#' @examples
#' \dontrun{
#' ## x is a list of sampled chains
#' calc_rhat(x)
#' }


calc_rhat <- function(x) {
  
  ## first, since the test cases for walkr passed
  ## we can safely assume here that x is of the correct format
  ## which we require it to be
  
  ## m is the number of chains
  
  m <- length(x)
  
  ## n is the number of points in each chain
  
  n <- dim(x[[1]])[2]
  
  ## params is the number of parameters (dimension of sampling space)
  ## that we have
  
  params <- dim(x[[1]])[1]
  rhats <- numeric()
  for (i in 1:params) {
    
    ## I'll document more later
    ## this is just some variance / mean calculation 
    ## across and within the different chains
    
    ## first, calculate the average value of each chain (miu bar_j)
    
    mu_each_chain <- as.numeric(lapply (x, function(y) { mean(y[i,])})   )
    
    ## the average value in the average value of chains 
    ## a constant which we need to calculate the variance across chains
    
    theta_2bar <- (1/m) * sum(mu_each_chain)
    
    ## see paper for details on B
    
    B <- (n / (m-1)) * sum (  (mu_each_chain - theta_2bar)^2 )
    
    ## see paper for details on W
    
    W <- sum(as.numeric(lapply(x, function(chain_matrix) {stats::var(chain_matrix[i, ])}))) / m
    
    ## calcualating Rhat squared
    
    R2 <- (W * (1 - 1/n) + B / n) / W
    
    ## returning the rhat of that parameter
    
    rhats <- c(rhats, sqrt(R2))
    
  }
  
  
  return(rhats)
  
  
  
}