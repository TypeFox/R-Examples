rgl.lines3d <- function (x, y, z, add = FALSE, ...) 
{
    n <- length(x)
    if (n != length(y) || n != length(z)) 
        stop("Lengths of x,  y, z do not match")
    if (!add) 
        rgl.clear()
    off <- c(1, 1)
    x <- kronecker(x, off)[-c(1, 2 * n)]
    y <- kronecker(y, off)[-c(1, 2 * n)]
    z <- kronecker(z, off)[-c(1, 2 * n)]

    rgl.lines(x,y,z,...)
}
