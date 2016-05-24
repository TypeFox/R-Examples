model.matrix.mars <-
function (object, x, which = object$selected.terms, full = FALSE, 
    ...) 
{
    if (missing(x)) 
        return(object$x)
    x <- as.matrix(x)
    dd <- dim(x)
    n <- dd[1]
    p <- dd[2]
    nterms <- length(which)
    dir <- object$factor
    cut <- object$cuts
    if (full) {
        bx <- matrix(0, nrow = n, ncol = object$lenb)
        bx[, 1] <- 1
    }
    else bx <- matrix(1, nrow = n, ncol = nterms)
    which <- which[-1]
    for (i in seq(along = which)) {
        j <- which[i]
        if (all(dir[j, ] == 0)) {
            stop("error in factor or which")
        }
        temp1 <- 1
        for (k in 1:p) {
            if (dir[j, k] != 0) {
                temp2 <- dir[j, k] * (x[, k] - cut[j, k])
                temp1 <- temp1 * temp2 * (temp2 > 0)
            }
        }
        if (full) 
            bx[, j] <- temp1
        else bx[, i + 1] <- temp1
    }
    bx
}

