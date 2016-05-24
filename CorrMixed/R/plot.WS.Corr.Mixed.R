plot.WS.Corr.Mixed <- plot.WS.Corr.Mixed.SAS <- 
  function (x, xlab, ylab, ylim, main, All.Individual=FALSE, ...){
  
if (missing(ylab)) {ylab="Reliability"}  
if (missing(ylim)) {ylim=c(-1, 1)}

  if (x$Model=="Model 1, Random intercept" | x$Model=="Model 1") {
    if (missing(xlab)) {xlab="Measurement moment"}
    if (missing(main)) {main="Model 1"}
    plot(y=x$CI.Upper, x=x$Time, xlab=xlab, ylab=ylab, ylim=ylim, col=0, main=main, ...)
    lines(y=x$R, x=x$Time)
    if (is.null(x$CI.Upper)==FALSE) {     
    lines(y=x$CI.Upper, x=x$Time, lty=2)
    lines(y=x$CI.Lower, x=x$Time, lty=2)}
  }

if (x$Model=="Model 2, Random intercept + serial corr (Gaussian)" | x$Model=="Model 2") {

  if (All.Individual==FALSE){
    if (missing(xlab)) {xlab="Time lag"}
    if (missing(main)) {main="Model 2"}
    plot(y=x$CI.Upper, x=x$Time, xlab=xlab, ylab=ylab, ylim=ylim, col=0, main=main, ...)
    lines(y=x$R, x=x$Time)
    if (is.null(x$CI.Upper)==FALSE) { 
    lines(y=x$CI.Upper, x=x$Time, lty=2)
    lines(y=x$CI.Lower, x=x$Time, lty=2)}
  }
  
  if (All.Individual==TRUE){
    if (missing(xlab)) {xlab="Measurement moment"}
    if (missing(main)) {main="Model 2"}
    plot(x=x$Time, y=1:max(length(x$R)), col=0, ylim=ylim, ylab=ylab, 
         xlab=xlab, main=main, ...)
    aantal <- 
      max(length(x$R))
    for (i in 1: aantal){
      line_hier <- x$R[1:(aantal-1)]
      lines(y=line_hier, x=x$Time[i : (length(line_hier)+(i-1))], col="grey")
    }
  }
}

if (x$Model=="Model 3, Random intercept, slope + serial corr (Gaussian)") {
  
    if (missing(xlab)) {xlab="Measurement moment"}
    if (missing(main)) {main="Model 3"}
    plot(x=x$Time, y=x$R[,1], 
         col=0, ylim=ylim, ylab=ylab, 
         xlab=xlab, main=main, ...)
    aantal <- dim(x$R)[2]
    for (i in 1: aantal){
      lines(y=x$R[,i][(i) : (aantal)], x=x$Time[(i) : (aantal)], col="grey")
      }
}  

if (x$Model=="Model 3") {
  
  if (missing(xlab)) {xlab="Measurement moment"}
  if (missing(main)) {main="Model 3"}
  plot(x=x$Time, y=x$R[,1], 
       col=0, ylim=ylim, ylab=ylab, 
       xlab=xlab, main=main, ...)
  aantal <- dim(x$R)[2]
  for (i in 1: aantal){
    lines(y=x$R[,i][(i) : (aantal)], x=x$Time[(i) : (aantal)], col="grey")
    
  }
    
}

}