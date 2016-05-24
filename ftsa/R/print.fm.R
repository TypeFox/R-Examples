print.fm = function (x, ...) 
{
    if (class(x)[1] == "fm"){ 
        cat("Functional model\n")
        cat("\nCall:", deparse(x$call, 200), "\n")
        mean <- is.element("mean", colnames(x$basis))
        level <- is.element("level", colnames(x$basis))
        order <- dim(x$T)[2]
        cat(paste("\nMain effects:", ifelse(mean, "Mean", "Mean"), ifelse(level, 
            "Level", "")))
        cat(paste("\nOrder of interaction:", order))
        cat(paste("\n   y:", x$y$yname, "\n   x:", x$y$xname, "\n"))
    }
    else {
         stop("object is not a functional model.")
    }
}
