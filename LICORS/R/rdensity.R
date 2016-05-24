#' @title Generate random sample from density() or wKDE
#'
#' @description 
#' This function draws random samples given \code{data} and a
#' \code{\link[stats]{density}} estimate (or just providing the correct
#' bandwidth \code{bw}).
#' 
#' @param n number of samples
#' @param fhat an object of class \code{'\link[stats]{density}'}
#' for bandwidth selection (if \code{bw} is not explicitly provided as argument)
#' @param weights vector of weights. Same length as \code{data}.
#' Default \code{weights=NULL} - in this case equal weights for each point
#' @param data underlying sample of \code{fhat}  
#' @param kernel kernel choice for \code{fhat}. Default: \code{kernel='Gaussian'}. 
#' See \code{\link[stats]{density}} for other options.
#' @param bw choice of bandwidth. Default: \code{bw=fhat$bw}. 
#' Again see \code{\link[stats]{density}} for other options.
#' @keywords distribution nonparametric
#' #@export
#' @examples
#' set.seed(1923)
#' xx = c(rnorm(100, mean = 2), runif(100))
#' aa = density(xx)
#' plot(aa)
#' xx_sample =  rdensity(n=1000, fhat = aa, data = xx)
#' lines(density(xx_sample) , col = 2)
#' 

rdensity <- function(n = 100, data = NULL, fhat = NULL, 
                     bw = fhat$bw, 
                     weights = NULL, 
                     kernel = "Gaussian") {
  nn <- n
  if (is.null(weights)) {
    weights <- rep(1 / length(data), length(data))
  }
  if (is.null(bw)) {
    stop("You must provide a bandwidth 'bw'.")
  }
  if (any(kernel == c("Gaussian", "gaussian"))){
    sample_data <- rnorm(nn, 
                         mean = sample(data, size = nn, 
                                       replace = TRUE, prob = weights), 
                         sd = bw)
  } else {
    stop("No other kernel than 'Gaussian' is implemented yet.")
  }
  return(sample_data)
} 
