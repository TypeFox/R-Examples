##
##  i n p o l y g o n . R  Polygon Interior
##


inpolygon <- function(x, y, xp, yp, boundary = FALSE) {
    stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y),
              is.numeric(xp), is.numeric(yp), length(xp) == length(yp))

    # xmin <- min(xp); xmax <- max(xp)  # exclude all points 'far' outside
    # ymin <- min(yp); ymax <- max(yp)

    n  <- length(x)
    np <- length(xp)

    # Polygon must be closed
    if (xp[1] != xp[np] || yp[1] != yp[np]) {
        xp <- c(xp, xp[1])
        yp <- c(yp, yp[1])
        np <- np + 1
    }

    inpoly <- rep(FALSE, n)
    onpoly <- rep(FALSE, n)

    j <- np
    for (i in 1:np) {
        dxp <- xp[j] - xp[i]
        dyp <- yp[j] - yp[i]
        dist <- dxp * (y - yp[i]) - (x - xp[i]) * dyp

        idx1 <- ( ((yp[i] <= y & y < yp[j])  |
                   (yp[j] <= y & y < yp[i])) &
                  (0 < dist * dyp) )
        inpoly[idx1] <- !inpoly[idx1]

        idx2 <- ( ((yp[i] <= y & y <= yp[j]) | (yp[j] <= y & y <= yp[i])) &
                  ((xp[i] <= x & x <= xp[j]) | (xp[j] <= x & x <= xp[i])) &
                  (0 == dist | !dxp) )
        onpoly[idx2] <- TRUE

        j <- i
    }

    if (boundary) inpoly[onpoly] <- TRUE
    else          inpoly[onpoly] <- FALSE

    return(inpoly)
}
