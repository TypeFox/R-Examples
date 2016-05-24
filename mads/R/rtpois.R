#' Ramdonly generates values from a zero-truncated Poisson distribution
#'
#' Generates values from a zero-truncated Poisson distribution with mean
#' equal to that specified. It uses a look up table to check which value of 
#' lambda will give values with the requested mean. 
#'  
#' @param N number of values to randomly generate
#' @param mean mean of the generated values
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @importFrom stats runif
#' @importFrom stats dpois
#' @importFrom stats qpois
#'
rtpois <-
function(N, mean=NA){
  truncated.poisson.table <- data.frame(mean = 1:20, lambda = c(5.6469669879783e-05, 1.5936343341849, 2.82144264293862, 3.9207052500513, 4.96512876739183, 5.98490244735758, 6.99357578449238, 7.9973210436761, 8.99889023979162, 9.99954616955456, 10.9998163134426, 11.999926276307, 12.9999591161807, 14, 15, 16, 17, 18, 19, 20))
  #find corresponding lambda value for desired mean
  lambda <- truncated.poisson.table$lambda[mean] #can make this just a vector rather than a dataframe
  lambda <- ifelse(is.na(lambda), mean, lambda)
  #generate quantiles from a uniform distribution between the probability of getting a zero, given lambda, and 1.
  #endpoints are altered to ensure p is between these values and not equal to as this generates 0's or Infs.
  p <- runif(N, dpois(0, lambda)+1e-10, 1-1e-10)
  #find the smallest integer x such that P(X <= x) >= p
  tpois <- qpois(p, lambda)
  return(tpois)
} 

