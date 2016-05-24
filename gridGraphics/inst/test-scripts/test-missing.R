
library(gridGraphics)

n <- 7

primtest2 <- function(nas, na) {
    angle <- seq(0, 2*pi, length=n+1)[-(n+1)]
    y <- 0.5 + 0.4*sin(angle)
    x <- 0.5 + 0.4*cos(angle)
    if (any(nas))
        text(x[nas], y[nas], paste("NA", (1:n)[nas], sep=""), col="gray")
    x[nas] <- na
    y[nas] <- na
    polygon(x, y, col="light gray", border=NA)
    lines(x, y)
}

celltest <- function(r, c, nas, na) {
    plot.new()
    primtest2(nas, na)
}

cellnas <- function(i) {
    temp <- rep(FALSE, n)
    temp[i] <- TRUE
    temp[n-3+i] <- TRUE
    temp
}

missing1 <- function() {
    par(mfrow=c(2, 2), mar=rep(0, 4), pty="s")
    celltest(1, 1, rep(FALSE, n), NA)
    celltest(1, 2, cellnas(1), NA)
    celltest(2, 1, cellnas(2), NA)
    celltest(2, 2, cellnas(3), NA)
}

plotdiff(expression(missing1()), "missing-1")

plotdiffResult()

