## last modified June 2002

conditdat <- function(mixdat, k, conditsamples) 
{
    xname <- dimnames(mixdat)[[2]][1:2]
    cdt <- matrix(0, nrow = nrow(mixdat), ncol = k)
    rowmat <- matrix(conditsamples, ncol = (k + 1), byrow = TRUE)
    cdt[as.vector(rowmat[, 1]), ] <- rowmat[, -1]
    dat <- as.data.frame(cbind(mixdat[, 1:2], cdt))
    if (is.null(xname)) 
        xname <- c("X", "count")
    dimnames(dat)[[2]] <- c(xname, paste("C", 1:k, sep = ""))
    dat
}
