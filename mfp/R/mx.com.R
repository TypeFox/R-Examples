mx.com <- function(x, y)
{
    if(is.matrix(x) & is.matrix(y)) {
        if(sum(dim(x) == dim(y)) != 2)
            stop("matrices have different dimensions")
    }
    x <- as.vector(x)
    y <- as.vector(y)
    if(length(x) != length(y))
        stop("vector lengths unequal")
    na.same <- is.na(x) == is.na(y) # NA status of x & y
    na.diff <- sum(!na.same)
    if(!na.diff) {
        xy.same <- x[!is.na(x)] == y[!is.na(y)] 
    # number status of x & y
        xy.diff <- sum(!xy.same)
    }
    else xy.diff <- 1
    return(!xy.diff)
}
