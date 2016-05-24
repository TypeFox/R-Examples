fds <- function (x, y, xname, yname) 
{
    if (missing(xname)) 
        xname <- deparse(substitute(x))
    if (missing(yname)) 
        yname <- deparse(substitute(y))
    y <- as.matrix(y)
    if (length(x) != nrow(y)) 
        stop("Dimensions do not match")
    return(structure(list(x = x, y = y, xname = xname, 
        yname = yname), class = "fds"))
}
