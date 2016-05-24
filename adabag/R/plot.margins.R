plot.margins <-
function(x, y=NULL, ...){
  if(!((class(x)=="margins")))
    stop("x class should be margins")
  if(!((class(y)=="margins")|is.null(y)))
    stop("y class should be margins or NULL")
  
  plot(sort(x$margins), (1:length(x$margins))/length(x$margins), 
       type="l", xlim=c(-1,1),main="Margin cumulative distribution graph", xlab="m", ylab="% observations", col="blue3", lwd=2)
  
  abline(v=0, col="red",lty=2, lwd=2)  
  if(!is.null(y)) {
    lines(sort(y$margins), (1:length(y$margins))/length(y$margins), type="l", cex = .5 ,col="green", lwd=2)
    legend("topleft", c("test","train"), col = c("blue", "green"), lty=1, lwd=2)
  }
}
