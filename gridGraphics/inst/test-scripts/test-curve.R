
library(gridGraphics)

curve1 <- function() {
    plot(qnorm) # default range c(0, 1) is appropriate here,
                # but end values are -/+Inf and so are omitted.
}

curve2 <- function() {
    plot(qlogis, main = "The Inverse Logit : qlogis()")
    abline(h = 0, v = 0:2/2, lty = 3, col = "gray")
}

curve3 <- function() {
    curve(sin, -2*pi, 2*pi, xname = "t")
}

curve4 <- function() {
    curve(tan, xname = "t", add = NA,
          main = "curve(tan)  --> same x-scale as previous plot")
}

curve5 <- function() {
    par(mfrow = c(2, 2))
    curve(x^3 - 3*x, -2, 2)
    curve(x^2 - 2, add = TRUE, col = "violet")

    ## simple and advanced versions, quite similar:
    plot(cos, -pi,  3*pi)
    curve(cos, xlim = c(-pi, 3*pi), n = 1001, col = "blue", add = TRUE)
}

chippy <- function(x) sin(cos(x)*exp(-x/2))

curve6 <- function() {
    curve(chippy, -8, 7, n = 2001)
}

curve7 <- function() {
    plot (chippy, -8, -5)
}

curve8 <- function() {
    par(mfrow=c(2, 2))
    for(ll in c("", "x", "y", "xy"))
        curve(log(1+x), 1, 100, log = ll, sub = paste0("log = '", ll, "'"))
}

plotdiff(expression(curve1()), "curve-1")
# Strange difference in positioning of "i" in "qlogis" in title
# There is some heavy rasterisation occurring, but not sure why it
# should be occurring differently for 'graphics' vs 'grid' version.
# Anyway, upping the resolution fixes things (on my system).
plotdiff(expression(curve2()), "curve-2", density=200)
plotdiff(expression(curve3()), "curve-3")
plotdiff(expression(curve4()), "curve-4")
plotdiff(expression(curve5()), "curve-5")
plotdiff(expression(curve6()), "curve-6")
plotdiff(expression(curve7()), "curve-7")
plotdiff(expression(curve8()), "curve-8")

plotdiffResult()
