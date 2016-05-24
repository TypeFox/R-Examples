#Create SSNPoint Class
setClass("SSNPoint", representation(
		network.point.coords = "data.frame",
		point.coords = "matrix",
		point.data = "data.frame",
		points.bbox = "matrix",
		proj4string = "CRS"
	),
	validity <- function(object) {
		if(length(object@network.point.coords[,1]) != length(object@point.data[,1])) {
			stop("Number of network point coordinate system is not equal to the number of spatial points")
		}
		#else
		TRUE
	}
)

#Create SSNPoints Class - holds a list of SSNPoint objects and their IDs
setClass("SSNPoints", 
	representation(
		SSNPoints = "list",
		ID = "character"
	),
)


#CREATE SpatialStreamNetwork CLASS
setClass("SpatialStreamNetwork", 
	representation("SpatialLinesDataFrame",
		network.line.coords = "data.frame",
		obspoints = "SSNPoints",
		predpoints = "SSNPoints",
		path = "character"
	),
	validity <- function(object) {
		if(length(object@network.line.coords[,1]) != length(object@data[,1]))
		{
			stop("Number of network line coordinate system is not equal to the number of spatial lines")
		}
		TRUE
	}
)

# create a default plotting method for SpatialStreamNetwork objects
setMethod("plot", signature(x = "SpatialStreamNetwork", y = "missing"),
  function(x, y, ...) plot.SpatialStreamNetwork(x,...))

# create a default summary method for SpatialStreamNetwork objects
setMethod("summary", signature(object = "SpatialStreamNetwork"),
  function(object, ...) show(object,...))


