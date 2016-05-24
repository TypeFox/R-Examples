"icllde" <-
function (I, h, f, m, n.iterations = 10, x1, xm, right.limit = 10000, kernel="gaussian") 
{
    if (missing(x1)) 
        x1 <- min(I) - 4 * h
    if (missing(xm)) 
        xm <- max(c(I[, 1], I[I[, 2] < right.limit, 2])) + 4 * 
            h
    if (missing(m)) 
        m <- length(f)
    x <- seq(x1, xm, length = m)
    xdiff <- x[2] - x[1]
    if (missing(f)) 
        f <- rep(1, m)/(m * xdiff)
    n <- dim(I)[1]
    left <- I[, 1]
    right <- I[, 2]
    numgrid <- m
    gridpts <- x
    ker <- which(c("gaussian", "epanechnikov", "biweight") %in% kernel)
    for (i in seq(1,n)) {
        continue <- any((I[i,2]>=gridpts)&(I[i,1]<=gridpts))
        if (!continue) break()
        }
    if (continue){
        f0 <- f
        niter <- n.iterations
        z <- .Fortran("icllde", as.integer(n), as.double(left), as.double(right), 
        as.integer(numgrid), as.double(gridpts), as.double(f0), 
        as.double(h), as.integer(niter), as.double(f), as.integer(ker),
        PACKAGE = "ICE")
        names(z) <- c("n", "left", "right", "numgrid", "x", "f0", 
        "h", "niter", "y", "ker")
        z.out <- list(x = z$x, y = z$y)
        class(z.out) <- "IC"
        z.out
        }        else {print("Error: At least one interval does not contain any gridpoints")}
}

