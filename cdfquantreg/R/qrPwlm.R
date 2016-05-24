
#' @title Probability Weighted L-moment Skewness and Kurtosis
#' @aliases qrPwlm
#' @description  Calculate the skew and kurtosis statistics based on probability weighted moments, via simulation method.
#' 
#' @param x The vector of values for the calculation of  Skewness and Kurtosis.
#' @param n The number of samples drawn in the simulation. The higher this value, the greater accuracy.  
#' @param mu vector of means.
#' @param sigma vector of standard deviations.
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#'
#' @export
#' 
#' @details This function computes the L-moment measures of skew and kurtosis, which may be computed via linear combinations of probability-weighted moments (Greenwood, Landwehr, Matalas and Wallis, 1979). 
#' 
#' @return The tau3(skew) and tau4(kurtosis) values of the L-moment.
#' 
#' @references Greenwood, J. A., Landwehr, J. M., Matalas, N. C., & Wallis, J. R. (1979). Probability weighted moments: definition and relation to parameters of several distributions expressable in inverse form. Water Resources Research, 15(5), 1049-1054.
#' 
#' @examples
#' qrPwlm(n = 1000, mu = 0.5, sigma = 1, fd = 't2', sd = 't2')
qrPwlm <- function(x, n = NULL, mu = NULL, sigma=NULL, fd=NULL, sd=NULL) {
  
  if(is.null(n)){
    n <- length(x)
  }
  # random sample
  if(missing(x)){
    x <- rq(n, mu, sigma, fd, sd)
    x <- sort(x)
  }

  # generate probablity weights
  betas <- rep(0, 4)
  for (k in 0:3) {
    i <- k + 1
    sum <- 0
    for (j in seq(1, n)) {
      sum <- sum + choose(j - 1, k) * x[j]
    }
    betas[i] <- sum/(n * choose(n - 1, k))
  }
  
  # generate 4 moments
  
  L2 <- 2 * betas[2] - betas[1]
  L3 <- 6 * betas[3] - 6 * betas[2] + betas[1]
  L4 <- 20 * betas[4] - 30 * betas[3] + 12 * betas[2] - betas[1]
  t3 <- L3/L2
  t4 <- L4/L2
  
  c(t3, t4,betas)
}
 
