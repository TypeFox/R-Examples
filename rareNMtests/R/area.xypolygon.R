area.xypolygon <-
function (polly) {
    xp <- polly[,1]
    yp <- polly[,2]
    nedges <- length(xp)
    yp <- yp - min(yp)
    nxt <- c(2:nedges, 1)
    dx <- xp[nxt] - xp
    ym <- (yp + yp[nxt])/2
    sum(dx * ym)
}
