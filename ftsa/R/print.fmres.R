print.fmres = function (x, ...) 
{
    if (class(x)[1] == "fmres") {
        cat("Residuals of a functional model\n")
        cat("\nCall:", deparse(x$call, 200), "\n")
        cat(paste("\nResiduals:"))
        cat(paste("\n   y:", x$yname, "\n   x:", x$xname, 
            "\n"))
    }
    else {
        stop("object is not a functional time series model")
    }
}