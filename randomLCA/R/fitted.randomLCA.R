`fitted.randomLCA` <-
function(object,...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
    fitted <- data.frame(object$patterns,Freq=object$freq,fitted=object$fitted)
    fitted
 }

