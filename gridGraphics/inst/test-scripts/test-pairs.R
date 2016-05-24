
library(gridGraphics)

pairs1 <- function() {
    pairs(iris[1:4], main = "Anderson's Iris Data -- 3 species",
          pch = 21, bg = c("red", "green3", "blue")[unclass(iris$Species)])
}

pairs2 <- function() {
    ## formula method
    pairs(~ Fertility + Education + Catholic, data = swiss,
          subset = Education < 20, main = "Swiss data, Education < 20")
}

pairs3 <- function() {
    pairs(USJudgeRatings)
}

pairs4 <- function() {
    ## show only lower triangle (and suppress labeling for whatever reason):
    pairs(USJudgeRatings, text.panel = NULL, upper.panel = NULL)
}

pairs5 <- function() {
    ## put histograms on the diagonal
    panel.hist <- function(x, ...) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
    }
    pairs(USJudgeRatings[1:5], panel = panel.smooth,
          cex = 1.5, pch = 24, bg = "light blue",
          diag.panel = panel.hist, cex.labels = 2, font.labels = 2)
}

pairs6 <- function() {
    ## put (absolute) correlations on the upper panels,
    ## with size proportional to the correlations.
    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste0(prefix, txt)
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
    }
    pairs(USJudgeRatings, lower.panel = panel.smooth, upper.panel = panel.cor)
}

pairs7 <- function() {
    pairs(iris[-5], log = "xy") # plot all variables on log scale
}

pairs8 <- function() {
    pairs(iris, log = 1:4, # log the first four
          main = "Lengths and Widths in [log]", line.main=1.5, oma=c(2,2,3,2))
}

plotdiff(expression(pairs1()), "pairs-1", width=10, height=10)
plotdiff(expression(pairs2()), "pairs-2")
plotdiff(expression(pairs3()), "pairs-3", width=15, height=15)
plotdiff(expression(pairs4()), "pairs-4", width=15, height=15)
plotdiff(expression(pairs5()), "pairs-5")
plotdiff(expression(pairs6()), "pairs-6", width=15, height=15, density=72)
plotdiff(expression(pairs7()), "pairs-7", width=10, height=10)
plotdiff(expression(pairs8()), "pairs-8", width=12, height=12)

plotdiffResult()
