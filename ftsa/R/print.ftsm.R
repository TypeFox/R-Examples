print.ftsm = function (x, ...) 
{
    if (class(x)[1] == "ftsm"){
        cat("Functional time series model\n")
        cat("\nCall:", deparse(x$call, 200), "\n")
        mean <- is.element("mean", colnames(x$basis))
        level <- is.element("level", colnames(x$basis))
        order <- ncol(x$basis) - mean - level
        cat(paste("\nMain effects:", ifelse(mean, "Mean", ""), ifelse(level, 
            "Level", "")))
        cat(paste("\nOrder of interaction:", order))
        cat(paste("\n   y:", x$y$yname, "\n   x:", x$y$xname, "\n"))
    }
    else {
         stop("object is not a functional time series model") 
    }       
}

