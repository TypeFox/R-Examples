`summary.cddews` <-
function (object, ...)
{
    ctmp <- class(object)
    if (is.null(ctmp))
        stop("cddews has no class")
    else if (ctmp != "cddews")
        stop("cddews is not of class cddews")
    cat("Locally stationary two-dimensional wavelet decomposition structure\n")
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    cat("Levels: ", object$nlevels, "\n")
    cat("dimension of original image was: ", object$datadim[1],
        "x", object$datadim[2], "pixels.\n")
    cat("Filter family used: ", object$family, "Filter index (N):", object$filter, "\n")
    cat("Structure adopted:", object$structure, "\n")
    cat("Date: ", object$date, "\n")
}

