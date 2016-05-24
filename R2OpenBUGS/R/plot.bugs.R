plot.bugs <- function (x, display.parallel = FALSE, ...){
    mar.old <- par("mar")
    pty.old <- par(pty = "m")
    mfrow.old <- par("mfrow")
        layout(matrix(c(1,2),1,2))
        
    bugs.plot.summary (x, ...)
    bugs.plot.inferences (x, display.parallel, ...)
    header <- ""
    if(!is.null(x$model.file))
        header <- paste(header, "Bugs model at \"", x$model.file, "\", ", sep="")
    if(!is.null(x$program))
        header <- paste(header, "fit using ", x$program, ", ", sep="")
    header <- paste(header, x$n.chains, " chains, each with ",
        x$n.iter, " iterations (first ", x$n.burnin, " discarded)", sep = "")
    mtext(header, outer = TRUE, line = -1, cex = 0.7)
    par(pty = pty.old[[1]], mar = mar.old, mfrow = mfrow.old)
}
