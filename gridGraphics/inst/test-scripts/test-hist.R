
library(gridGraphics)

require(utils) # for str
require(stats)

hist1 <- function() {
    par(mfrow = c(2, 2))
    hist(islands)
    utils::str(hist(islands, col = "gray", labels = TRUE))

    hist(sqrt(islands), breaks = 12, col = "lightblue", border = "pink")
    ##-- For non-equidistant breaks, counts should NOT be graphed unscaled:
    r <- hist(sqrt(islands), breaks = c(4*0:5, 10*3:5, 70, 100, 140),
              col = "blue1")
    text(r$mids, r$density, r$counts, adj = c(.5, -.5), col = "blue3")
    lines(r, lty = 3, border = "purple") # -> lines.histogram(*)
}

hist2 <- function() {
    hist(islands, breaks = c(12,20,36,80,200,1000,17000), freq = TRUE,
         main = "WRONG histogram") # and warning
}

hist3 <- function() {
    set.seed(14)
    x <- rchisq(100, df = 4)
    par(mfrow = 2:1, mgp = c(1.5, 0.6, 0), mar = .1 + c(3,3:1))
    ## Comparing data with a model distribution should be done with qqplot()!
    qqplot(x, qchisq(ppoints(x), df = 4)); abline(0, 1, col = 2, lty = 2)

    ## if you really insist on using hist() ... :
    hist(x, freq = FALSE, ylim = c(0, 0.2))
    curve(dchisq(x, df = 4), col = 2, lty = 2, lwd = 2, add = TRUE)
}

plotdiff(expression(hist1()), "hist-1", width=10, height=10)
plotdiff(expression(hist2()), "hist-2")
plotdiff(expression(hist3()), "hist-3")

plotdiffResult()
