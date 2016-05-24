gx.ngr.ltdl <-
function (xmat, vars = NULL, coded = -9999) 
{
    if (!(is.matrix(xmat) | is.data.frame(xmat))) 
        stop(paste("  ", deparse(substitute(xmat)),
            "is not a matrix or data frame"))
    if.df <- FALSE
    if (is.data.frame(xmat)) {
        xsav <- xmat
        ind.num <- sapply(xmat, is.numeric)
        xmat <- as.matrix(xmat[, ind.num])
        if.df = TRUE
    }
#
    if(is.null(vars)) {
        xname <- dimnames(xmat)[[2]]
        nvars <- dim(xmat)[2]
        vars <- seq(1:nvars)
    }
    else {
        nvars <- length(vars)
        xname <- character(nvars)
        varnums <- integer(nvars)
        for (i in 1:nvars) {
            ii <- vars[i]
            if (is.numeric(vars[i])) xname[i] <- dimnames(xmat)[[2]][ii]
            else xname[i] <- vars[i]
        }
    }
    cat("  Variable        N     NA     <DL (-ve)", coded, "\n\n")
#
    for (i in 1:nvars) {
        ii <- vars[i]
        x <- xmat[, ii]
        n <- length(x)
        ncoded <- length(x[!is.na(x) & x == coded])
        x[x == coded] <- 0
        nna <- length(x[is.na(x)])
        nneg <- length(x[!is.na(x) & x < 0])
        cat(" ", xname[i], "    \t", n, "\t", nna, "\t", nneg, "\t   ",
            ncoded, "\n")
    }
    invisible()
}
