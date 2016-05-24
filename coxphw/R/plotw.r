plotw <- function
(
  x,                      # object of class coxphw
  rank=FALSE,
  log=FALSE,
  legendxy=NULL,
  ...                     
)
{
  w.matrix<-x$w.matrix[order(x$w.matrix[,1]),]
  if(rank) {
    time<-rank(w.matrix[,1])
    label<-"ranked time"
  }
  else {
    time<-w.matrix[,1]
    label<-"time"
  }
  if(log) {
    weights<-log(w.matrix[,2:4])
    wlabel<-"log of weight"
  }
  else {
    weights<-w.matrix[,2:4]
    wlabel<-"weight"
  }
  
  plot(time, weights[,3],type="l",lty=1,ylim=c(min(weights),max(weights)), xlab=label, ylab=wlabel, 
       lwd=2, las=1, ...)   
  lines(time, weights[,1],lty=2)                                                      
  lines(time, weights[,2],lty=3)
  if (is.null(legendxy)) { legendxy <- c(min(time),0.95*max(weights)) }
  
  if (x$template %in% "AHR") {
  legend(x=legendxy[1], y=legendxy[2], title="weights", 
         legend=c("survival","censoring", "combined (normalized)"), lty=c(2:3, 1), lwd=c(1,1,2), bty="n")
  } else
  if (x$template %in% "ARE") {
    legend(x=legendxy[1], y=legendxy[2], title="weights", 
           legend=c("survival","censoring", "combined (normalized)"), lty=c(2:3, 1), lwd=c(1,1,2), bty="n")
  } else {
    legend(x=legendxy[1], y=legendxy[2], lty=1:3, lwd=c(2,1,1), bty="n",
           legend=c("combined (normalized)", "raw weight","censoring weight"))
  }
}
