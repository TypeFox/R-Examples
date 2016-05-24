factDesign <- function (dimension, levels) 
{
    if (length(levels) == 1) {
        levels <- rep(levels, dimension)
    }
    else if (length(levels) != dimension) {
        stop("argument 'levels' must be of length 1 or 'dimension'")
    }
    x <- matrix(0, prod(levels), dimension)
    for (i in 1:dimension) {
        xi <- 1
        for (j in (1:dimension)[(1:dimension) < i]) {
            xi <- xi %x% rep(1, levels[j])
        }
        xi <- xi %x% seq(0, 1, length.out = levels[i])
        for (j in (1:dimension)[(1:dimension) > i]) {
            xi <- xi %x% rep(1, levels[j])
        }
        x[, i] <- xi
    }
    n <- dim(x)[1]
    return(list(design = x, n = n, dimension = dimension, levels = levels))
}