## last modified June 2002

mixgroup <- function(x, breaks = NULL, xname = NULL, k = NULL, usecondit = FALSE) 
{
    if (is.data.frame(x) | is.matrix(x)) {
        xx <- x[, 1]
        if (is.null(xname)) 
            xname <- dimnames(x)[[2]][1]
    }
    else xx <- x
    if (is.null(breaks)) {
        xh <- hist(xx, plot = FALSE)
    }
    else {
        xh <- hist(xx, breaks = breaks, plot = FALSE)
    }
    xhdf <- data.frame(X = c(xh$breaks[c(-1, -length(xh$breaks))], 
        Inf), count = xh$counts)
    if (!is.null(xname)) 
        dimnames(xhdf)[[2]][1] <- xname
    if (usecondit & (is.data.frame(x) | is.matrix(x)) & ncol(as.matrix(x)) > 
        1) {
        itab <- table(cut(xx, xh$breaks), x[, 2])
        attr(itab, "class") <- "matrix"
        kk <- as.numeric(dimnames(itab)[[2]][ncol(itab)])
        if (is.null(k)) 
            k <- kk
        else if (k < kk) {
            warning(paste(kk, "groups found, argument k =", k, 
                "ignored"))
            k <- kk
        }
        xmat <- matrix(0, nrow = nrow(itab), ncol = k)
        dna <- as.numeric(dimnames(itab)[[2]])
        xmat[, dna] <- itab
        xhdf <- cbind(xhdf, xmat)
        dimnames(xhdf)[[2]][3:(k + 2)] <- c(paste("C", 1:k, sep = ""))
    }
    class(xhdf) <- c("mixdata", "data.frame")
    xhdf
}
