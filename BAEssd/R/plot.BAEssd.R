plot.BAEssd <-
function(x,y="TE",alpha.line=TRUE,type="l",
    xlab="Sample Size (n)",ylab=NULL,main=NULL,...){
  
  ### Error checking
  if(!(y %in% colnames(x$history))){
    stop("Argument y is not a column of history.")
  }
  
  ### Determine the y-axis label if not provided
  if(is.null(ylab)){
    ylab <- switch(y,
        TE="Total Error",
        TWE="Total Weighted Error",
        AE1="Average Bayes Type-I Error",
        AE2="Average Bayes Type-II Error")
  }
  
  ### Create plot
  plot(x$history[,y]~x$history[,"n"],type=type,
      xlab=xlab,ylab=ylab,main=main,...)
  
  if(alpha.line) abline(h=attr(x$n,"alpha"),lty=2)
}
