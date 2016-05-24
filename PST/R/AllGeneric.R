## 

setGeneric(
      name="cmine",
	def=function(object, ...)
	standardGeneric("cmine")
)

setGeneric(name="cplot",
	def=function(object, context, ...)
	standardGeneric("cplot")
)

setGeneric("cprob", 
	def=function(object, L, ...) 
	standardGeneric("cprob")
)

setGeneric("gain", 
	def=function(object, ...) 
	standardGeneric("gain")
)

setGeneric("generate", 
	def=function(object, ...) 
	standardGeneric("generate")
)

setGeneric(
      name="impute",
	def=function(object, data, ...)
	standardGeneric("impute")
)

setGeneric(
      name="nodenames",
	def=function(object, ...)
	standardGeneric("nodenames")
)

setGeneric(
      name="pdist",
	def=function(x, y, ...)
	standardGeneric("pdist")
)

setGeneric(
      name="pmine",
	def=function(object, data, ...)
	standardGeneric("pmine")
)

setGeneric(
      name="ppplot",
	def=function(object, path, ...)
	standardGeneric("ppplot")
)

setGeneric(
      name="pqplot",
	def=function(object, data, ...)
	standardGeneric("pqplot")
)


setGeneric(
      name="prune",
	def=function(object, ...)
	standardGeneric("prune")
)

setGeneric(
      name="pstree",
	def=function(object, group, L, ...)
	standardGeneric("pstree")
)

setGeneric("query", 
	def=function(object, context, state, output="prob", exact=FALSE, ...)
	standardGeneric("query")
)

setGeneric(
      name="subtree",
	def=function(object, ...)
	standardGeneric("subtree")
)

setGeneric(name="tune",
	def=function(object, ...)
	standardGeneric("tune")
)



