setGeneric(name="apcluster",
           def=function(s, x, ...) standardGeneric("apcluster"))

setGeneric(name="apclusterL",
           def=function(s, x, ...) standardGeneric("apclusterL"))

setGeneric(name="apclusterK",
           def=function(s, x, ...) standardGeneric("apclusterK"))

setGeneric(name="aggExCluster",
           def=function(s, x, ...) standardGeneric("aggExCluster"))

setGeneric(name="heatmap",
           def=function(x, y, ...) standardGeneric("heatmap"))

setGeneric(name="similarity",
           def=function(x, ...) standardGeneric("similarity"))

setGeneric(name="preferenceRange",
           def=function(s, ...) standardGeneric("preferenceRange"))

setGeneric(name="as.SparseSimilarityMatrix",
           def=function(s, ...) standardGeneric("as.SparseSimilarityMatrix"))

setGeneric(name="as.DenseSimilarityMatrix",
           def=function(s, ...) standardGeneric("as.DenseSimilarityMatrix"))
