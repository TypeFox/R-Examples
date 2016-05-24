## ----plot----------------------------------------------------------------
library(RcppFaddeeva)
x <- seq(-1000, 1000)
x0 <- 200
l <- Lorentz(x, x0, 30)
g <- Gauss(x, x0, 100)
N <- length(x)
c <- convolve(Gauss(x, 0, 100), 
              rev(Lorentz(x, x0, 30)), type="o")[seq(N/2, length=N)]
v <- Voigt(x, x0, 100, 30)
matplot(x, cbind(v, l, g, c), t="l", lty=c(1,2,2,1), xlab="x", ylab="")
legend("topleft", legend = c("Voigt", "Lorentz", "Gauss", "Convolution"), bty="n",
      lty=c(1,2,2,1), col=1:4)

## ----integrals-----------------------------------------------------------
integrate(Lorentz, -Inf, Inf, x0=200, gamma=100)
integrate(Gauss, -Inf, Inf, x0=200, sigma=50)
integrate(Voigt, -Inf, Inf, x0=200, sigma=50, gamma=100)

## ----complex-------------------------------------------------------------
x <- seq(-1000, 1000)
x0 <- 200
v <- Voigt(x, x0, 100, 30, real = FALSE)
matplot(x, cbind(Re(v), Im(v)), t="l", lty=c(1,2), xlab="x", ylab="", col=1)
legend("topleft", legend = c("Imaginary part", "Real part"), bty="n",
      lty=c(1,2), col=1)

