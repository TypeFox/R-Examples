alr <-
function (xx, j = NULL, ifclose = FALSE, ifwarn = TRUE) 
{
    if (!is.matrix(xx)) 
        stop("  ", deparse(substitute(xx)), " is not a Matrix")
    temp.x <- remove.na(xx)
    x <- temp.x$x
    p <- temp.x$m
    if (ifwarn) 
        cat("  ** Are the data all in the same measurement units? **\n")
    if (is.null(j)) 
        stop("  ** The divisor must be specified **")
    if (j > p) 
        stop("j cannot be >", p)
    if (ifclose) 
        x <- 100 * sweep(x, 1, rowSums(x), "/")
    x <- log(x)
    x <- sweep(x, 1, x[, j], "-")
    x <- as.matrix(x[, -j])
    return(x = x)
}
