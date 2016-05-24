gdt <- function(y) {
  # Author: Francisco Parra Rodr?guez 
  # http://econometria.wordpress.com/2013/08/21/estimation-of-time-varying-regression-coefficients/ 
  a <- matrix(y,nrow=1)
  n <- length(y)
  A <- MW(n)
  t(A)%*%t(a)
}