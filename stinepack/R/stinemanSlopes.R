"stinemanSlopes" <-
function (x, y, scale = FALSE) 
{
    m <- length(x)
    m1 <- m - 1
    if (m == 2) {
        yp <- rep(diff(y)/diff(x), 2)
    }
    else {
        if (scale) {
            sx <- diff(range(x))
            sy <- diff(range(y))
            if (sy <= 0) 
                sy <- 1
            x <- x/sx
            y <- y/sy
        }
        dx <- diff(x)
        dy <- diff(y)
        yp <- rep(NA, m)
        dx2dy2p <- dx[-1]^2 + dy[-1]^2
        dx2dy2m <- dx[-m1]^2 + dy[-m1]^2
        yp[2:m1] <- (dy[-m1] * dx2dy2p + dy[-1] * dx2dy2m)/(dx[-m1] * 
            dx2dy2p + dx[-1] * dx2dy2m)
        s <- dy[1]/dx[1]
        if ((s >= 0 && s >= yp[2]) || (s <= 0 && s <= yp[2])) 
            yp[1] <- 2 * s - yp[2]
        else yp[1] <- s + abs(s) * (s - yp[2])/(abs(s) + abs(s - 
            yp[2]))
        s <- dy[m1]/dx[m1]
        if ((s >= 0 && s >= yp[m1]) || (s <= 0 && s <= yp[m1])) 
            yp[m] <- 2 * s - yp[m1]
        else yp[m] <- s + abs(s) * (s - yp[m1])/(abs(s) + abs(s - 
            yp[m1]))
        if (scale) 
            yp <- yp * sy/sx
    }
    yp
}
