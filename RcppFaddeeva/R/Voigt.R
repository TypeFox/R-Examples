##' @description Voigt distribution
##' @title The Voigt function, corresponding to the convolution of a lorentzian and a gaussian distribution
##' @param x numeric vector
##' @param x0 scalar, peak position
##' @param sigma parameter of the gaussian
##' @param gamma parameter of the lorentzian
##' @param real logical, return only the real part of the complex Faddeeva
##' @param ... passed to Faddeeva_w
##' @return numeric or complex vector
##' @author baptiste Auguie
##' @describeIn Voigt Voigt lineshape function
##' @family helper_function
##' @examples
##' ## should integrate to 1 in all cases
##' integrate(Lorentz, -Inf, Inf, x0=200, gamma=100)
##' integrate(Gauss, -Inf, Inf, x0=200, sigma=50)
##' integrate(Voigt, -Inf, Inf, x0=200, sigma=50, gamma=100)
##' 
##' ## visual comparison
##' x <- seq(-1000, 1000)
##' x0 <- 200
##' l <- Lorentz(x, x0, 30)
##' g <- Gauss(x, x0, 100)
##' N <- length(x)
##' c <- convolve(Gauss(x, 0, 100), 
##'               rev(Lorentz(x, x0, 30)), type="o")[seq(N/2, length=N)]
##' v <- Voigt(x, x0, 100, 30)
##' matplot(x, cbind(v, l, g, c), t="l", lty=c(1,2,2,1), xlab="x", ylab="")
##' legend("topleft", legend = c("Voigt", "Lorentz", "Gauss", "Convolution"), bty="n",
##'        lty=c(1,2,2,1), col=1:4)
##' @export
Voigt <- function(x, x0, sigma, gamma, real = TRUE, ...){
  
  z <- (x - x0 + gamma*1i) / (sigma * sqrt(2))
  w <- Faddeeva_w(z, ...)
  if(real) return(Re(w) / (sigma * sqrt(2*pi))) else
    w / (sigma * sqrt(2*pi))
  
}

##' @description Lorentzian distribution
##' @title Lorentz
##' @describeIn Voigt Lorentzian lineshape function
##' @inheritParams Voigt 
##' @export
##' @family helper_function
Lorentz <- function(x, x0, gamma){
  gamma / (pi*((x -x0)^2 + gamma^2))
}

##' @description Gaussian distribution
##' @title Gauss 
##' @describeIn Voigt Gaussian lineshape function
##' @inheritParams Voigt
##' @export
##' @family helper_function
Gauss <- function(x, x0, sigma){
  dnorm(x, x0, sd = sigma)
}



