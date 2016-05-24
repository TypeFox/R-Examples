
library(gridGraphics)

mtext1 <- function() {
    plot(1:10, (-4:5)^2, main = "Parabola Points", xlab = "xlab")
    mtext("10 of them")
    for(s in 1:4)
        mtext(paste("mtext(..., line= -1, {side, col, font} = ", s,
                    ", cex = ", (1+s)/2, ")"), line = -1,
              side = s, col = s, font = s, cex = (1+s)/2)
    mtext("mtext(..., line= -2)", line = -2)
    mtext("mtext(..., line= -2, adj = 0)", line = -2, adj = 0)
}

mtext2 <- function() {
    ##--- log axis :
    plot(1:10, exp(1:10), log = "y", main = "log =\"y\"", xlab = "xlab")
    for(s in 1:4) mtext(paste("mtext(...,side=", s ,")"), side = s)
}

plotdiff(expression(mtext1()), "mtext-1", width=10, height=10)
plotdiff(expression(mtext2()), "mtext-2", width=10, height=10)

plotdiffResult()
