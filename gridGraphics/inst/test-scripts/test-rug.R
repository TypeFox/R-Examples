
library(gridGraphics)

require(stats)  # both 'density' and its default method

rug1 <- function() {
    with(faithful, {
        plot(density(eruptions, bw = 0.15))
        rug(eruptions)
        rug(jitter(eruptions, amount = 0.01), side = 3, col = "light blue")
    })
}

plotdiff(expression(rug1()), "rug-1", antialias=FALSE)

plotdiffResult()
