loess.sd <- function (x, y = NULL, nsigma = 1, ...)
{
# Loess mean smooth +- k*sd
# based on a similar function available on Arc
#
# Arguments:
# x, y = data 
# k = multiplier for standard deviation
# ... = further parameters passed to loess()
#
# Example:
# data(cars)
# plot(cars, main = "lowess.sd(cars)")
# lines(l <- loess.sd(cars))
# lines(l$x, l$upper, lty=2)
# lines(l$x, l$lower, lty=2)
#
# Reference:
# Weisberg, S. (2005) Applied Linear Regression, 3rd ed., Wiley, New York, 
#   pp. 275-278.

  xy <- xy.coords(x, y)
  x <- xy$x
  x0 <- sort(x)
  y <- xy$y
  nsigma <- as.numeric(nsigma)
  mod <- loess(y ~ x, ...)
  yfit <- predict(mod, data.frame(x = x0))
  r  <- residuals(mod)
  modr <- loess(I(r^2) ~ x, ...)
  sd <- sqrt(pmax(0, predict(modr, data.frame(x = x0))))
  list(model = mod, x = x0, y = yfit, sd = sd, 
       upper = yfit + nsigma*sd, lower = yfit - nsigma*sd)
}

panel.loess <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "red", span = 2/3, degree = 2, nsigma = 1, ...) 
{
# Panel smoothing using loess with variability bands at nsigmas*sd
# See 'loess.sd'
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
       { l <- loess.sd(x[ok], y[ok], nsigma = nsigma, 
                       span = span, degree = degree, ...)
         lines(l$x, l$y, col = col.smooth, ...)
         lines(l$x, l$upper, lty = 2, col = col.smooth, ...)
         lines(l$x, l$lower, lty = 2, col = col.smooth, ...) 
       }
}
