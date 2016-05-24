setMethod("length", signature(x="APResult"), function(x) length(x@exemplars))

setMethod("length", signature(x="AggExResult"), function(x) x@maxNoClusters)

setMethod("length", signature(x="ExClust"), function(x) length(x@exemplars))
