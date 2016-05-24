##
##  f m i n v i z . R  Function Minimum and Line Test
##


fminviz <- function(fn, x0, nlines = 2*length(x0),
                    npoints = 51, scaled = 1.0) {
    n <- length(x0)
    R <- matrix(runif(n * nlines, min = -1, max = 1),
                nrow = n, ncol = nlines)
    for (i in 1:nlines)
        R[, i] <- R[, i] / sqrt(sum(R[, i] * R[, i]))

    N <- npoints
    x <- seq(-1, 1, length.out = N) * scaled
    Y <- matrix(0, nrow = N, ncol = nlines)
    for (j in 1:nlines) {
        for (i in 1:N) {
            Y[i, j] <- fn(x0 + x[i] * R[, j])
        }
    }
    plot(scaled * c(-1, 1), c(min(Y), max(Y)), type = "n",
         xlab = "Distance (from x0)", ylab = "Function values",
         main = "Function Minimum Test")
    grid()

    y0 <- fn(x0)
    for (j in 1:nlines) {
        if (all(Y[, j] >= y0)) clr <- "blue"
        else                   clr <- "red"
        lines(x, Y[, j], col = clr)
    }
    invisible(NULL)
}


flineviz <- function(fn, x1, x2, npoints = 51, scaled = 0.1) {
    n <- length(x1)
    if (length(x2) != n)
        stop("")

    N <- npoints
	ii <- seq(-0.1, 1.1, length.out = N)
	ll <- sum((x2-x1)^2)
	xs <- ys <- numeric(N)

	for (i in 1:N) {
	    xs[i] <- ii[i]*ll
	    ys[i] <- fn(x1 + ii[i]*(x2 - x1))
	}

	plot(xs, ys, type="l", lwd=2, col="gray",
	     xlab = "Distance (between x1 and x2)", ylab = "Function values",
	     main = "Function Line Test")
	points(c(0, ll), c(fn(x1), fn(x2)), col = "red")
	grid()

	invisible(NULL)
}



