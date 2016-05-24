BIC.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (missing(...)) {
        if (object$method == "REML") {
            object$fit.stats$REML[4]
        }
        else {
            object$fit.stats$ML[4]
        }
    }
    else {
        if (object$method == "REML") {
            out <- sapply(list(object, ...), function(x) x$fit.stats$REML[4])
        }
        else {
            out <- sapply(list(object, ...), function(x) x$fit.stats$ML[4])
        }
        dfs <- sapply(list(object, ...), function(x) x$parms)
        out <- data.frame(df = dfs, BIC = out)
        Call <- match.call()
        rownames(out) <- as.character(Call[-1L])
        return(out)
    }
}
