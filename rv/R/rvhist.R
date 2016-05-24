# ========================================================================
# rvhist  -  plot histograms of the simulations of the random components
# ========================================================================
#

rvhist <- function (x, ...)  {
    if (!is.null(dim(x))) 
        par(mfcol = dim(x))
    mfcol <- par("mfcol")
    n <- prod(mfcol)
    n <- min(n, length(x))
    a <- list(...)
    if (is.null(a$freq) && is.null(a$prob)) {
      a$freq <- FALSE
    }
    if (is.null(a$breaks)) {
      a$breaks <- "fd"
    }
    make.main <- is.null(a$main)
    make.xlab <- is.null(a$xlab)
    lab <- deparse(substitute(x))
    x.names <- paste(lab, .dimindex(x), sep = "")
    for (i in 1:n) {
        a$x <- sims(x[i])
        if (make.main || make.xlab) {
            this.name <- x.names[i]
            if (make.xlab) {
                a$xlab <- this.name
            }
            if (make.main) {
                a$main <- paste("Histogram of", this.name)
            }
        }
        do.call("hist", a)
    }
}


