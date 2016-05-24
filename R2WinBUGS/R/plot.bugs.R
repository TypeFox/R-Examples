plot.bugs <- function (x, display.parallel = FALSE, ...){
    mar.old <- par("mar")
    pty.old <- par(pty = "m")
    mfrow.old <- par("mfrow")
    if (is.R())
        layout(matrix(c(1,2),1,2))
    else
        par(mfrow = c(1,2))
        
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
    if (is.R())  par(pty = pty.old[[1]], mar = mar.old, mfrow = mfrow.old)
    else  invisible(par(pty = pty.old[[1]], mar = mar.old, mfrow = mfrow.old))
}

if (!is.R()) {

strwidth <-function(s, units = c("user", "inches", "figure"), cex = NULL) {
   s<-as.character(s)
   if (!missing(cex)) {
     oldcex <- par(cex=cex)
     on.exit(par(oldcex))
   }
   units <- match.arg(units)
   if (units == "user") {
      nchar(s) * par("cxy")[1]
   } else if (units == "inches") {
      nchar(s) * par("cin")[1]
   } else if (units == "figure") {
      nchar(s) * par("cin")[1] / par("fin")[1]
   }
}

strheight <- function(s, units = "user", cex = NULL) {
   s<-as.character(s)
   if (!missing(cex)) {
     oldcex <- par(cex=cex)
     on.exit(par(oldcex))
   }
   units <- match.arg(units)
   if (units == "user") {
      par("cxy")[2]
   } else if (units == "inches") {
      par("cin")[2]
   } else if (units == "figure") {
      par("cin")[2] / par("fin")[2]
   }
}

} #ends if (!is.R())
