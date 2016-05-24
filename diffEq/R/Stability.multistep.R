## =============================================================================
## Generates Stability functions for Multistep and Runge-Kutta methods
## =============================================================================
# ------------------------------------------------------------------------------
# Utility function
# given a vector 'x' and a set of polynomial coefficients 'coef', 
# calculates the polynomial c[1] + c[2]*x + c[3]*x^2 + c[4]*x^3 ... 
# for all elements in x
# ------------------------------------------------------------------------------

polynom <- function (x, coef = 1) {
  p  <- 0:(length (coef)-1)
  xx <- matrix(ncol = length(x), nrow = length(coef), 
               data = x, byrow = TRUE)
  colSums(xx^p*coef)
}      

## -----------------------------------------------------------------------------
## Stability plotting function for Multistep methods
## -----------------------------------------------------------------------------

multistepBnd <- function (alpha, 
                          beta) {

  theta <- seq(0, 2*pi, by = 0.01)
  z     <- exp(1i*theta)

  polynom(z,alpha)/polynom(z,beta)
  
}

stability.multistep<- function (alpha, 
                          beta, 
                          add = FALSE, 
                          fill = NULL,
                          ...) {
  nu <- multistepBnd(alpha, beta)
  if (! add ) 
    plot(nu, type = "l", xlab = "Re(z)", ylab = "Im(z)", ...)
  else
    lines(nu, ...)

  if (! is.null(fill)) {
    polygon(nu, col = fill, border = "black")
  }
  abline(h = 0)
  abline(v = 0)

}