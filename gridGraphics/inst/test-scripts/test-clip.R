
library(gridGraphics)

clip1 <- function() {
    set.seed(1)
    x <- rnorm(1000)
    hist(x, xlim = c(-4,4))
    usr <- par("usr")
    clip(usr[1], -2, usr[3], usr[4])
    hist(x, col = 'red', add = TRUE)
    clip(2, usr[2], usr[3], usr[4])
    hist(x, col = 'blue', add = TRUE)
}

plotdiff(expression(clip1()), "clip-1")

plotdiffResult()
