AIC.rma <-
function (object, ..., k = 2, correct = FALSE) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (missing(...)) {
        if (object$method == "REML") {
            ifelse(correct, object$fit.stats$REML[5], object$fit.stats$REML[3])
        }
        else {
            ifelse(correct, object$fit.stats$ML[5], object$fit.stats$ML[3])
        }
    }
    else {
        if (object$method == "REML") {
            out <- sapply(list(object, ...), function(x) ifelse(correct, 
                x$fit.stats$REML[5], x$fit.stats$REML[3]))
        }
        else {
            out <- sapply(list(object, ...), function(x) ifelse(correct, 
                x$fit.stats$ML[5], x$fit.stats$ML[3]))
        }
        dfs <- sapply(list(object, ...), function(x) x$parms)
        out <- data.frame(df = dfs, AIC = out)
        if (correct) 
            names(out)[2] <- "AICc"
        Call <- match.call()
        Call$k <- NULL
        Call$correct <- NULL
        rownames(out) <- as.character(Call[-1L])
        return(out)
    }
}
