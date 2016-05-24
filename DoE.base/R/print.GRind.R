print.GRind <- function(x, quote=FALSE, ...){
   GR.i <- x$GR.i
   GR.i <- as.character(GR.i)
   if (!is.null(GR.i)) GR.i[x$GR.i>=floor(x$GRs[2])+1] <- paste(">=",floor(x$GRs[2])+1,sep="")
   names(GR.i) <- names(x$GR.i)
   x$GR.i <- GR.i
   print.default(x, quote=quote, ...)
}