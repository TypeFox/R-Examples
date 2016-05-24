`print.cddews` <-
function (x, ...)
{
    if (IsEarly(x)) {
        ConvertMessage()
        stop()
    }
    cat("Class 'cddews' : corrected directional dependent wavelet spectrum:\n")
    cat("       ~~~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    if (x$correct == TRUE) {
        cat("The spectrum of this image was corrected (IP matrix).\n")
        cat("$S is a large array of data \n")
    }
    else {
        cat("This spectrum is UNCORRECTED.\n")
        cat("$spec is a large array of data of dimension", dim(x$spec),
            ".\n")
    }
    if (x$smooth == TRUE) {
        cat("The spectra have been smoothed to obtain consistency. \n")
    }
    else cat("The spectra were NOT smoothed!!! n")
    cat("\nCreated on :", x$date, "\n")
}


