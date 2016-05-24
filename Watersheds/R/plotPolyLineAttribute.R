
# S4 generic method plot for
# plotting attributes of a SpatialLinesDataFrame or
# SpatialPolygonsDataFrame objects
# e.g. plots river or watershed OBJECTID

setGeneric("plot.PolyLineAttribute", function(x,y,dist,cex,...) standardGeneric("plot.PolyLineAttribute"))

setMethod(
	"plot.PolyLineAttribute", signature=c("SpatialPolygonsDataFrame", "character", "numeric", "numeric"), 
	function(x,y,dist,cex,...){
		midpoint = gCentroid(x, byid=TRUE)
		aa = midpoint@coords
		labels =  slot(x, "data") 
		labels = as.matrix(labels[y])
		text(aa+dist, col="blue", label=labels, cex=cex)	
		}	
)

setMethod(
	"plot.PolyLineAttribute", signature=c("SpatialLinesDataFrame", "character", "numeric", "numeric"), 
	function(x,y,dist,cex,...){
		midpoint = gCentroid(x, byid=TRUE)
		aa = midpoint@coords
		labels =  slot(x, "data") 
		labels = as.matrix(labels[y])
		text(aa+dist, col="blue", label=labels, cex=cex)	
		}	
)

