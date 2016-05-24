plot.CI <- function(x,...){
  
  # histogram of posterior samples
  hist(x,100,freq=FALSE,
       main=paste("posterior distribution for ",deparse(substitute(x)),sep=" "),
       xlab=deparse(substitute(x)),...)
  
  # indicate 95% credible interval, mean, and median
  abline(v=quantile(x,c(.025,.975)),lwd=2,col="red")
  abline(v=mean(x),lwd=2,col="green")
  abline(v=median(x),lwd=2,col="blue")
  
  # legend
  legend("topright",
         legend=c(paste(c("lower bound 95% CI","upper bound 95% CI"),
                        round(quantile(x,c(.025,.975)),3),sep=": "),
                  paste("mean",round(mean(x),3),sep=": "),
                  paste("median",round(median(x),3),sep=": ")),
         col=c("red","red","green","blue"),lty=1,lwd=2)
}