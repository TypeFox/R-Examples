"[.cnv" <-
function(obj, i) {
  old.attr <- attributes(obj)
  probs <- attr(obj,"probabilities")[i,,drop=FALSE]
  meanRatio<- if (!is.null(attr(obj,"meanRatio"))) .subset(attr(obj,"meanRatio"),i) else NULL
  batches<- if (!is.null(attr(obj,"batches"))) .subset(attr(obj,"batches"),i) else NULL
  obj <- .subset(obj, i)
  attributes(obj)<-old.attr
  if (!is.null(batches)){
    bb <- sort(unique(batches))
    attr(obj, "means") <- attr(obj, "means")[bb,,drop=FALSE]
    attr(obj, "sds") <- attr(obj, "sds")[bb,,drop=FALSE]
    attr(obj, "pi") <- attr(obj, "pi")[bb,,drop=FALSE]  
  }
  attr(obj,"meanRatio")<-meanRatio
  attr(obj,"batches")<-batches
  attr(obj,"probabilities")<-probs
  obj
}
