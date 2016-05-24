multhist <- function (x, beside=TRUE, freq=NULL, probability=!freq,
 plot.it=TRUE, ...) {
  ## sort out histogram arguments
  hist.args <- formals(hist.default)
  args <- list(...)
  hargs <- names(args)[names(args) %in% names(hist.args)]
  hist.args[hargs] <- args[hargs]
  hist.args$plot<-FALSE
  ## sort out barplot arguments
  barplot.args <- formals(barplot.default)
  bargs <- names(args)[names(args) %in% names(barplot.args)]
  barplot.args[bargs] <- args[bargs]
  barplot.args$beside <- beside
  ## prevent warnings
  barplot.args$"..." <- barplot.args$inside <- NULL
  allhist <- hist(unlist(x),hist.args$breaks,plot=FALSE)
  if (!"names.arg" %in% bargs) {
    barplot.args$names.arg <- signif(allhist$mids, 2)
  }
  if (is.null(freq)) {
   freq<-if(!missing(probability)) !as.logical(probability)
    else TRUE
  }
  if(freq) comp<-"counts"
  else comp<-"density"
  combhist <- t(sapply(x,
   function(z) hist(z,breaks=allhist$breaks,plot=FALSE)[[comp]]))
  if(plot.it) do.call("barplot", c(list(combhist),barplot.args))
  invisible(list(allhist,combhist))
}
