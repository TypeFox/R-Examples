rkoch <-
function (n, side = 3, niter = 5) 
{
    vert <- koch(side, niter)
    square <- c(min(vert[, 1]), max(vert[, 1]), min(vert[, 2]), 
        max(vert[, 2]))
    m = 0
    xy <- matrix(nrow = n, ncol = 2)
    while (m < n) {
        x <- runif(1, square[1], square[2])
        y <- runif(1, square[3], square[4])
        if (in.polygon(x, y, vert[, 1], vert[, 2])) {
            xy[m + 1, ] <- c(x, y)
            m <- m + 1
        }
    }
    return(xy)
}
