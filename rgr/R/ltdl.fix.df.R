ltdl.fix.df <-
function (x, negs2na = FALSE, zero2na = FALSE, coded = NA) 
{
    if (!(is.matrix(x) | is.data.frame(x))) 
        stop(paste("  ", deparse(substitute(x)), "is not a matrix or data frame"))
    if.df <- FALSE
    n.f <- 0
    if (is.data.frame(x)) {
        xsav <- x
        ind.num <- sapply(x, is.numeric)
        x <- as.matrix(x[, ind.num])
        if.df <- TRUE
        n.f <- length(ind.num[ind.num == FALSE])
    }
    p <- length(x[1, ])
    n <- length(x[, 1])
    nna <- sum(is.na(x))
    cat("  n =", n, "by p =", p, "matrix checked,", nna, "NA(s) present\n ", 
        n.f, "factor variable(s) present")
    ncoded <- 0
    if (!is.na(coded)) {
        x[x == coded] <- NA
        ncoded <- sum(is.na(x)) - nna
        cat("\n ", ncoded, "value(s) coded", coded, "set to NA")
    }
    if (zero2na) {
        x[abs(x) < 10^-5] <- NA
        nzero <- sum(is.na(x)) - nna - ncoded
        cat("\n ", nzero, "zero (abs(x) < 10^-5) value(s) set to NA")
    }
    nfix <- length(x[!is.na(x) & x < 0])
    if (negs2na) {
        x[!is.na(x) & x < 0] <- NA
        cat("\n ", nfix, "-ve value(s) set to NA\n")
    }
    else {
        x[!is.na(x) & x < 0] <- abs(x[!is.na(x) & x < 0])/2
        cat("\n ", nfix, "-ve value(s) set to +ve half the negative value\n")
    }
    x <- as.data.frame(x)
    if (if.df) {
        xsav[, ind.num] <- x
        x <- xsav
    }
}
