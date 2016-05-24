summary.checks<-function(object, ...) {
  list(Means=object@means,items=colMeans(object@tab,na.rm=TRUE))
}
