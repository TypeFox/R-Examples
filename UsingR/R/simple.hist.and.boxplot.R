##' Plot both the histogram and the boxplot  to show relationship easily.
##'
##' @param x data
##' @param ... passed to hist
##' @return NULL
##'
##' @export
"simple.hist.and.boxplot" <-
  function (x,...) {
    op<-par(no.readonly=TRUE)
    on.exit(par(op))
    layout(matrix(c(1,2),2,1),heights=c(3,1))
    par(mai=c(1,1,1,1)/2)
    hist(x,xlab=FALSE,col=gray(0.95),yaxt='n',...)
    rug(x)
    boxplot(x,horizontal=TRUE)
  }
