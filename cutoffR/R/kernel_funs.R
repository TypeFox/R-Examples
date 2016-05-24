#' The Uniform Kernel
#' 
#' @param x function arguments
#' @export
#' @examples
#' 
#' curve(UnifK)
#' plot(UnifK, -2, 2)
UnifK <- function(x) ifelse(abs(x) <=1, 0.5, 0)

#' The Cosine Kernel
#' 
#' @param x function arguments
#' @export
#' @examples
#' 
#' curve(CosK)
#' plot(CosK, -2, 2)
CosK <- function(x) ifelse(abs(x) <= 1, pi/4*cos(pi*x/2), 0)

#' The Epanechnikov Kernel
#' 
#' @param x function arguments
#' @export
#' @examples
#' 
#' curve(EpanK)
#' plot(EpanK, -2, 2)
EpanK <- function(x) ifelse(abs(x) <= 1, 3/4*(1-x^2), 0)


#' The Gaussian Kernel
#' 
#' @param x function arguments
#' @export
#' @examples
#' 
#' curve(GaussK)
#' plot(GaussK, -2, 2)
GaussK <- function(x) 1/sqrt(2*pi)*exp(-1/2*x^2)
