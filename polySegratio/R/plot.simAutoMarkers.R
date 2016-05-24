`plot.simAutoMarkers` <-
function(x, main=deparse(substitute(x)),
                                xlab="Segregation ratio", ...) {

  ## description: plot segregation ratios of

  ## Arguments:
  ## x: object of class simAutoMarkers
  ## main: title for plot
  ## xlab: label for x axis

  main # hmm is this for Splus
##  plot(hist(rowMeans(x$markers,na.rm = TRUE)),
##       main=paste("Segregation ratios for",main), xlab=xlab, ...)
  plot(x$seg.ratio, main=main,...)
}

