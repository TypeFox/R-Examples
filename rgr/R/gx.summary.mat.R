gx.summary.mat <-
function (xmat, vars, banner = deparse(substitute(xmat)), log = FALSE) 
{
    if (!(is.matrix(xmat) | is.data.frame(xmat))) 
        stop(paste("  ", deparse(substitute(xmat)), "is not a matrix or data frame"))
    if.df <- FALSE
    if (is.data.frame(xmat)) {
        xsav <- xmat
        ind.num <- sapply(xmat, is.numeric)
        xmat <- as.matrix(xmat[, ind.num])
        if.df <- TRUE
    }
    nvars <- length(vars)
    if (log) 
        cat("  Data log10 transformed: SD, CV% and SE in log10 units\n")
    cat("  Summary Stats for", banner, "\n\n\t", "N NAs - Min Q1 M Q2 Max - MAD IQR_SD - Mean SD CV% - SE 95% CI on Mean\n")
    for (i in 1:nvars) {
        ii <- vars[i]
        if (is.numeric(vars[i])){ 
            if (if.df) xname <- dimnames(xmat)[[2]][ii]
            else {xname <- as.character(ii)}
           }
        else {
            xname <- vars[i]
        }
        table <- gx.summary(xmat[, ii], log = log, iftell = FALSE)
        cat("\n  ", xname, ":\t ", table[1], " ", table[2], " - ", 
            table[3], " ", table[4], " ", table[5], " ", table[6], 
            " ", table[7], " - ", table[8], " ", table[9], " - ", 
            table[10], " ", table[11], " ", table[12], "% - ", 
            table[13], " ", table[14], " <-> ", table[15], "\n", 
            sep = "")
    }
    cat("\n")
    invisible()
}
