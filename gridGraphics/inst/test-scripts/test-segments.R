
library(gridGraphics)

segments1 <- function() {
    set.seed(1)
    x <- stats::runif(12); y <- stats::rnorm(12)
    i <- order(x, y); x <- x[i]; y <- y[i]
    plot(x, y, main = "arrows(.) and segments(.)")
    ## draw arrows from point to point :
    s <- seq(length(x)-1)  # one shorter than data
    arrows(x[s], y[s], x[s+1], y[s+1], col= 1:3)
    s <- s[-length(s)]
    segments(x[s], y[s], x[s+2], y[s+2], col= 'pink')
}

plotdiff(expression(segments1()), "segments-1")

plotdiffResult()
