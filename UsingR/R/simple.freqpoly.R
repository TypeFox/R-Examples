##' make a simple freqpoly and histogram
##'
##' @param x data
##' @param ... passed to hist
##' @return NULL
##'
##' @export
"simple.freqpoly" <-
function (x,...) {

tmp<-hist(x,probability=FALSE,...)
lines(c(min(tmp$breaks),tmp$mids,max(tmp$breaks)),c(0,tmp$counts,0),type="l")
}
