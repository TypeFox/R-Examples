GCV <- function (lambda,y,x,d) {
  #***************************************************************************
  # This function calculates the generalized cross validation (GCV) criterion
  #
  #   Input:  Lamda - value of lambda to calculate the GCV criterion at
  #           x - the design matrix for the penalized spline regression (PSR)
  #           d - the diagonal matrix used for smoothing the PSR
  #
  #   References:
  #         
  #   Hutchinson, Michael F., and F. R. De Hoog. "Smoothing noisy data with 
  #   spline functions." Numerische Mathematik 47.1 (1985): 99-106. 
  #
  #***************************************************************************
  s <- x%*%solve(t(x)%*%x+lambda^2*d)%*%t(x)
  mod <- lm(y~s%*%y)
  RSS <- sum(resid(mod)^2)
  GCV <- RSS/(1 - 1/length(y)*sum(diag(s)))^2
  return(GCV)
}