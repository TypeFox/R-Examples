plot.iemg <- function(x, type = "l", timeunits = c("seconds", "samples"), reset.lty = 2, 
    add = FALSE, ...) {
    object <- x
    args <- list(...)
    namesargs <- names(args)
    timeunits <- match.arg(timeunits)
    if ((timeunits == "seconds") & (object$samplingrate <= 0)) 
        warning("Sampling rate must be provided to determine the total duration of the iEMG")
    
    x <- 1:length(object$values)
    if ((timeunits == "seconds") & (object$samplingrate > 0)) {
        x <- x/object$samplingrate
        xlab = "time (seconds)"
    } else {
        xlab = "time (samples)"
    }
    tx <- x[object$reset.points]
    x[object$reset.points] <- NA
    
    if (!("ylab" %in% namesargs)) {
        if (object$data.name != "") 
            ylab <- object$data.name else ylab <- "iEMG"
        if (object$units != "") 
            ylab <- paste(ylab, "(", object$units, ")", sep = "")
        if (add) 
            lines(x, object$values, type = type, ylab = ylab, ...) else plot(x, object$values, type = type, ylab = ylab, xlab = xlab, ...)
    } else {
        if (add) 
            lines(x, object$values, type = type, ...) else plot(x, object$values, type = type, xlab = xlab, ...)
    }
    x[!(object$reset.points + 1)] <- NA
    x[object$reset.points] <- tx
    lines(head(x, length(object$values)), object$values, lty = reset.lty)
} 
