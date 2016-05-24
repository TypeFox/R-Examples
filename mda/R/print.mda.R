print.mda <-
function (x, ...) 
{
    NextMethod("print")
    if (!is.null(x$deviance)) 
        cat("\nDeviance:", format(round(x$deviance, 3)), "\n")
    invisible(x)
}

