
library(gridGraphics)

(wwt <- hist(women$weight, nclass = 7, plot = FALSE))

plot.histogram1 <- function() {
    plot(wwt, labels = TRUE) # default main & xlab using wwt$xname
}

plot.histogram2 <- function() {
    plot(wwt, border = "dark blue", col = "light blue",
         main = "Histogram of 15 women's weights", xlab = "weight [pounds]")
    ## Fake "lines" example, using non-default labels:
    w2 <- wwt; w2$counts <- w2$counts - 1
    lines(w2, col = "Midnight Blue", labels = ifelse(w2$counts, "> 1", "1"))
}

plotdiff(expression(plot.histogram1()), "plot.histogram-1")
plotdiff(expression(plot.histogram2()), "plot.histogram-2")

plotdiffResult()
