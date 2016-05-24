
# S4 generic method plot for
# plotting attributes of a SpatialPointDataFrame 
# e.g. plots node elevations or rivers Strahler

setGeneric("plot.PointAttribute", function(x,y,dist,cex,...) standardGeneric("plot.PointAttribute"))

setMethod(
	"plot.PointAttribute", signature=c("SpatialPointsDataFrame", "character", "numeric", "numeric"), 
	function(x,y,dist,cex,...){
		#x1=slot(x, "x"); x1
		#y1=slot(y, "y"); y1
		#dist1=slot(dist, "dist"); dist1
		#cex1=slot(cex, "cex")
		a = slot(x, "coords"); a
		labels =  slot(x, "data") 
		labels = as.matrix(labels[y])
		text(x=a[,1]+dist, y=a[,2]+dist, col="red", label=labels, cex=cex)
	}
)
