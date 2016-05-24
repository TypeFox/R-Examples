"plot.tune.ahazpen" <- function(x, df=TRUE, ...)
{
    ## Purpose: Plot 'tune.ahazpen' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'tune.ahazpen' object
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

  if(x$tune$type=="CV") {
    ylabname<-"Cross-validation score"
    yrng<-range(c(x$tuneup,x$tunelo))
  } else {
     ylabname<-"PBIC"
    yrng<-range(x$tunem)
  }

  # Check for updated arguments
  plot.args<-list(x=x$lambda,y=x$tunem,ylim=yrng,ylab=ylabname,log = "x", xlab = expression(lambda),pch = 20, col = 2)
  new.args<-list(...)
    if(length(new.args))
      plot.args[names(new.args)]<-new.args
   do.call("plot",plot.args)

  # Mark minumum lambda
  abline(v = x$lambda.min,lty = 2)

  # Degrees of freedom
  if(df) {
    vv<-(c(1,diff(x$df))!=0)
    axis(side = 3, at = x$lambda[vv], labels = paste(x$df[vv]), tick = TRUE, line = 0)
    mtext("# nonzero coefficients", side = 3, line = 2, adj = .5, cex = 0)
  }
  # Error bars
  if(x$tune$type=="CV")
    error.bars(x$lambda, x$tuneup, x$tunelo, width = 0.001, col = "darkgrey")
}
