"plot.cumahaz"<-function(x, ...)
  {
    ## Purpose: plot cumulative hazard function
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object: 'cumahaz' object
    ##   ...   : additional arguments to plot function
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    mint <- min(x$times)
    maxt <- x$times[order(x$times,decreasing=TRUE)[2]]#max(which(x$event != 0))]
    plot.args <- list(x=stepfun(x$times, x$cumhaz),main="",xlab="Time",
                    ylab="Cumulative baseline hazard",
                    do.points = FALSE, xlim = c(mint, maxt))
    new.args <- list(...)
    if(length(new.args))
      plot.args[names(new.args)] <- new.args
    do.call("plot.stepfun",plot.args)
  }
