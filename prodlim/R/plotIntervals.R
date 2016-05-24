plotIntervals <- function(object,rightCensored=FALSE,xlim,ylim,ylab,xlab,...){
  stopifnot(match("Hist",class(object)) && attr(object,"cens.type")=="intervalCensored")
  x <- object[order(object[,"L"]),]
  L <- x[,"L"]
  R <- x[,"R"]
  if (rightCensored==FALSE)
    x <- x[!is.na(R)&!is.infinite(R),]
  if (missing(ylim))
    ylim <- c(0,NROW(x)+1)
  if (missing(xlim))
    xlim <- c(0,max(R[!is.na(R)&!is.infinite(R)]))
  if (missing(xlab))
    xlab <- "Time"
  if (missing(ylab))
    ylab <- "Observed intervals"
  plot(0,0,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  nix <- lapply(1:NROW(x),function(f){
    x <- unlist(x[f,c("L","R"),drop=TRUE])
    x[is.infinite(x)] <- max(R[!is.na(R)&!is.infinite(R)])
    segments(x0=x[1],y0=f,x1=x[2],y1=f,lwd=2)})
  invisible(x)
}
