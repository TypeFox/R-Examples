##
##  r e c t i n t . R  Rectangular Intersections
##


rectint <- function(x, y) {
    stopifnot(is.numeric(x), is.numeric(y))
    if (is.vector(x) && length(x) == 4 &&
        is.vector(y) && length(y) == 4) {

        if (any(c(x[3], x[4], y[3], y[4]) < 0))
            stop("All widths and heights must be greater than 0.")

        if (x[1]+x[3] <= y[1] || y[1]+y[3] <= x[1] ||
            x[2]+x[4] <= y[2] || y[2]+y[4] <= x[2]) {
            return(0)

        } else {
            if (x[1] > y[1]) {
                tmp <- x; x <- y; y <- tmp
            }
            z1 <- y[1]
            z2 <- max(x[2], y[2])
            z3 <- min(x[1]+x[3], y[1]+y[3])
            z4 <- min(x[2]+x[4], y[2]+y[4])
            area <- (z3-z1) * (z4-z2)
            return(area)
        }

    } else if (is.matrix(x) && ncol(x) == 4 &&
               is.matrix(y) && ncol(y) == 4) {

        nx <- nrow(x); ny <- nrow(y)
        R <- matrix(NA, nrow = nx, ncol = ny)

        for (i in 1:nx) {
            for (j in 1:ny) {
                R[i, j] <- rectint(x[i, ], y[j, ])
            }
        }

        return(R)  
      
    } else {
        stop("All lengths and no. of matrix columns must be equal to 4.")
    }
}
