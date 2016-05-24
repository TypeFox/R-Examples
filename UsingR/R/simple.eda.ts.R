##' functin to make simple set of graphs for a times series data set
##'
##' @param x data
##' @param lag what lag
##' @return NULL
##'
##' @export
"simple.eda.ts" <-
  function(x,lag=1) {
    op <- par(no.readonly = TRUE);on.exit(par(op))
    par(mfrow=c(1,3)) 
    plot(x,main="Run sequence plot")
    plot(x,x[c((1+lag):length(x),1:lag)],main=paste("lag plot, lag =",lag))
#    library(ts)
    acf(x)
  }
