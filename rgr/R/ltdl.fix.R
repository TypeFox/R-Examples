ltdl.fix <-
function (x, negs2na = FALSE, zero2na = FALSE, coded = NA) 
{
    if (!is.numeric(x)) {
        cat("\nThe function argument must be numeric\n")
        return()
    }
    n <- length(x)
    nna <- sum(is.na(x))
    cat(" ", n, "records checked,", nna, "NA(s) present")
    ncoded <- 0
    if (!is.na(coded)) {
        x[x == coded] <- NA
        ncoded <- sum(is.na(x)) - nna
        cat("\n ", ncoded, "value(s) coded", coded, "set to NA")
    }
    if (zero2na) {
        x[abs(x) < 10^-5] <- NA
        nzero <- sum(is.na(x)) - nna - ncoded
        cat("\n ", nzero, "zero (abs(x) < 10^-5) record(s) set to NA")
    }
    nfix <- length(x[!is.na(x) & x < 0])
    if (negs2na) {
        x[!is.na(x) & x < 0] <- NA
        cat("\n ", nfix, "-ve record(s) set to NA\n")
    }
    else {
        x[!is.na(x) & x < 0] <- abs(x[!is.na(x) & x < 0])/2
        cat("\n ", nfix, "-ve record(s) set to +ve half the negative value\n")
    }
    invisible(x)
}
