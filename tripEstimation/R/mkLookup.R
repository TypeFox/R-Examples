"mkLookup" <- function (x, by.segment = TRUE)
{

    if (any(is.na(x$z))) stop("NAs in grid data")
    if (!by.segment & !is.logical(x$z)) stop("grid data must be a logical matrix")

    csize <- c(diff(x$x[1:2]), diff(x$y[1:2]))
    dimXY <- dim(x$z)
    binArray <- FALSE
    if (length(dimXY) == 3 & by.segment) {
        bsegs <- (1:(dimXY[3] * 31)%/%31) * prod(dimXY[1:2])
        dimXY <- dimXY[1:2]
        binArray <- TRUE
    }
    function(xy, segment = 1:nrow(xy)) {
        xs <- xy[, 1]
        ys <- xy[, 2]

        i <- round((1/diff(x$x[1:2]))*(xs - x$x[1])+1)
        j <- round((1/diff(x$y[1:2]))*(ys - x$y[1])+1)

        f <- vector(mode(x$z), length(xs))
        k <- (i > 0 & j > 0 & i <= dimXY[1] & j <= dimXY[2])
        n <- nrow(xy)
        if (any(k)) {
            if (binArray) {
                f[k] <- bits(x$z[((j[k] - 1) * dim(x$z)[1] +
                  i[k]) + bsegs[1:n][k]], (segment[k] - 1)%%31)
                f == 1
            }
            else {
                f[k] <- x$z[cbind(i[k], j[k])]
                f == 1
            }
        }
        else FALSE
    }
}
