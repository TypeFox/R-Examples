`cts` <-
function(x,y){
    if (!is.ts(x))
        stop("error: x must be ts object")
    if (nargs() != 2)
        stop("error: exactly two arguments required")
    if (is.ts(y))
        if(!all(end(x)+c(1,0) == start(y)))
            warning("y coerced to vector and tsp attribute ignored")
    f <- tsp(x)[3]
    y <- as.vector(y)
    ny <- length(y)
    if (ny == 0) return(x)
    endz <- end(x)+c(floor(ny/f),ny%%f)
    if (endz[2] > f)
        endz <- endz + c(1, -f)
    z <- window(x, start=start(x), end=endz, extend=TRUE)
    nz <- length(z)
    z[(nz-ny+1):nz]<-y
    z
}
