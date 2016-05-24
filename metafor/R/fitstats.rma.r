fitstats.rma <-
function (object, ..., REML) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (missing(REML)) {
        if (object$method == "REML") {
            REML <- TRUE
        }
        else {
            REML <- FALSE
        }
    }
    if (missing(...)) {
        if (REML) {
            out <- cbind(object$fit.stats$REML)
            colnames(out) <- "REML"
        }
        else {
            out <- cbind(object$fit.stats$ML)
            colnames(out) <- "ML"
        }
    }
    else {
        if (REML) {
            out <- sapply(list(object, ...), function(x) x$fit.stats$REML)
        }
        else {
            out <- sapply(list(object, ...), function(x) x$fit.stats$ML)
        }
        out <- data.frame(out)
        Call <- match.call()
        Call$REML <- NULL
        names(out) <- as.character(Call[-1L])
    }
    rownames(out) <- c("logLik:", "deviance:", "AIC:", "BIC:", 
        "AICc:")
    return(out)
}
