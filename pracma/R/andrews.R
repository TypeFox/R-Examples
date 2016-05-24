##
##  a n d r e w s . R  Andrews Curves
##


andrewsplot <- function(A, f, style = "pol", scaled = FALSE, npts = 101) {
    stopifnot(is.numeric(A), is.matrix(A))
    if (is.factor(f)) f <- as.integer(f)
    if (!is.integer(f))
        stop("Argument 'f' must be a factor or an integer vector.")

    n <- nrow(A); m <- ncol(A)
    if (m < 2 || n < 2)
        stop("Matrix 'A' must have at least two rows and columns.")

    if (scaled) A <- scaled(A)

    xpts <- seq(0, 2*pi, length = npts)
    Y <- matrix(NA, nrow = n, ncol = npts)

    # Compute the Andrews function
    for (i in 1:n) {
        xs <- A[i, ]
        ypts <- c()
        for (p in xpts) {
            y <- xs[1]
            for (j in 2:m) {
                if (j %% 2 == 1) { y <- y + xs[j]*sin((j %/% 2)*p) }
                else             { y <- y + xs[j]*cos((j %/% 2)*p) }
            }
            ypts <- c(ypts, y)
        }
        Y[i, ] <- ypts
    }

    if (style == "cart") {
        # plot in cartesian coordinates
        ymin <- min(Y)
        ymax <- max(Y)
        plot(0, 0, type="n", xlim=c(0, 2*pi), ylim=c(ymin, ymax),
             main="Andrews' Curves", xlab="", ylab="")
        for (i in sample(1:n, n)) {
            lines(xpts, Y[i, ], col = f[i])
        }

    } else if (style == "pol") {
        ymax <- max(abs(Y))
        polar(0, ymax, type="n", main = "Andrews' Curves", bxcol = "white")
        for (i in sample(1:n, n)) {
            polar(xpts, Y[i, ], col=f[i], add=TRUE)
        }

    } else
        stop("Argument 'style' can only be 'cart' or 'pol'.")

    invisible()
}
