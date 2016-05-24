ellipse <-
function(hlaxa = 1, hlaxb = 1, theta = 0, xc = 0, 
        yc = 0, npoints = 100, col) {
        xp <- NULL
        yp <- NULL
        for (i in 0:npoints) {
            a <- (2 * pi * i)/npoints
            x <- hlaxa * cos(a)
            y <- hlaxb * sin(a)
            if (theta != 0) {
                alpha <- angle(x, y)
                rad <- sqrt(x^2 + y^2)
                x <- rad * cos(alpha + theta)
                y <- rad * sin(alpha + theta)
            }
            xp <- c(xp, x + xc)
            yp <- c(yp, y + yc)
        }
        lines(xp, yp, type = "l", col = col)
        invisible()
    }
