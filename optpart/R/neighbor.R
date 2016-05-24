neighbor <- function (x, all = FALSE) 
{
    if (inherits(x, "pam")) {
        tmp <- silhouette(x)
        if (!all) 
            tmp <- tmp[tmp[, 3] < 0, ]
        out <- table(tmp[, 1], tmp[, 2])
    } else if (inherits(x, "partana")) {
        y <- x$clustering
        if (!all) {
            z <- apply(x$ptc,1,which.max)
            out <- table(y[y!=z],z[y!=z])
        } else {
            x <- x$ptc
            for (i in 1:nrow(x)) {
                x[i,y[i]] <- NA
            }
            z <- apply(x,1,which.max)
            out <- table(y,z) 
        }
    } else {
        stop("The first argument must be of class pam or partana")
    }
    out
}
