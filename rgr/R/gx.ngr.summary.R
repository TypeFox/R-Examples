gx.ngr.summary <-
function (xmat, vars = NULL, banner = deparse(substitute(xmat)), file = NULL) 
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
    if(is.null(file)) folder <- getwd()
    else folder <- file
    if (if.df) prefix <- "dataframe"
    else prefix <- banner
    filename = paste(folder, "/", prefix, "_NGR_summary.csv", sep = "")
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
    rownames <- c("N","NAs","Mean","Std. Dev.","Skewness","CV %","Geom Mean",
        "Median","MAD","Robust CV %","Minimum","1st %ile","2nd %ile",
        "5th %ile","10th %ile","20th %ile","25th %ile","30th %ile",
        "40th %ile","50th %ile","60th %ile","70th %ile","75th %ile",
        "80th %ile","90th %ile","95th %ile","98th %ile","99th %ile","Maximum")
    cat("  NGR Summary Stats for:\n  ", xname, "\n  from:", banner, "\n\n ",
        "NGR summary statistics will be saved in:\n  ", filename, "\n\n")
#
    table <- matrix(0, nvars, 29)
    for (i in 1:nvars) {
        ii <- vars[i]
        temp <- gx.ngr.stats(xmat[, ii])
        for (j in 1:29) table[i, j] <- temp[j]
    }
    for.csv <- matrix(0, 30, nvars+1)
    for.csv[1,1] <- " "
    for (j in 1:nvars) {
        for.csv[1, j+1] <- xname[j]
        for (i in 1:29) {
            for.csv[i+1, 1] <- rownames[i]
            for.csv[i+1, j+1] <- table[j, i]
        }
    }
    write.csv(for.csv, file = filename, row.names = FALSE)
    invisible()
}
