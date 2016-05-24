plot.eaemg <- function(x, ci.lty = "dashed", ...) {
    object <- x
    args <- list(...)
    namesargs <- names(args)
    if ("lty" %in% namesargs) {
        lty <- args$lty
    } else {
        lty <- "solid"
    }
    if (!("xlab" %in% namesargs)) {
        args <- c(args, xlab = "Percent cycle")
    }
    if (!("ylab" %in% namesargs)) {
        args <- c(args, ylab = "Prediction intervals")
    }
    
    if (is.null(object$intervals) || !is.matrix(object$intervals)) 
        stop("invalid 'envelope averaged emg' structure")
    n <- dim(object$intervals)[1]
    Percent.cycle <- matrix(tail(seq(0, 100, by = 100/n), n), nrow = n, ncol = 3)
    do.call("matplot", c(list(x = Percent.cycle, y = object$intervals, type = "l", 
        col = "black", lty = c(ci.lty, lty, ci.lty)), args))
} 
