subset.cnv <-
function(x,subset, ...) {
  old.attr<-attributes(x)
  probs<-subset.matrix(attr(x,"probabilities"),subset)
  if (!is.null(attr(x,"meanRatio"))) meanRatio<-subset.default(attr(x,"meanRatio"),subset)
  x<-subset.default(x,subset)
  attributes(x)<-old.attr
  if (!is.null(attr(x,"meanRatio"))) attr(x,"meanRatio")<-meanRatio
  attr(x,"probabilities")<-probs
  x
}

