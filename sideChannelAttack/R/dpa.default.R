dpa1.default <-
function (x, y, ...) 
{
    if (!is.matrix(x)) {
        stop("'x' has to be a matrix")
        return(-1)
    }
    if (!is.factor(y)) {
        stop("'y' has to be a factor")
        return(-1)
    }
    if (dim(x)[1] != length(y)) {
        stop("the number of output has to be the same as the number of input")
        return(-1)
    }
    if (length(which(y==0)) > 0) {
        stop("'y' has to contain values strictly greater than 0")
        return(-1)
    }
    MoyenneDeLa <- list()
    for (cle in sort(unique(y))) {
        MoyenneDeLa[[cle]] = apply(matrix(x[which(y == cle), ], ncol = dim(x)[2]), MARGIN = 2, FUN = mean)
    }
    res = list(mean = MoyenneDeLa)
    class(res) <- "dpa1"
    return(res)
}