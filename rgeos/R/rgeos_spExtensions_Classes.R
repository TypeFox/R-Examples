#########################
###
### SpatialRings
###
#########################

setClass("Ring", 
    representation(coords = "matrix",ID = "character"),
    validity = function(object) {
        if (any(is.na(object@coords)))
            stop("coords cannot contain missing values")
        if (ncol(object@coords) != 2)
            stop("coords should have 2 columns")
        if (!all(object@coords[1,] == object@coords[nrow(object@coords),]))
            stop("invalid ring, start and end coordinates must be equal")
        return(TRUE)
    }
)


setClass("SpatialRings",
    representation("Spatial", rings = "list"),
    
    prototype = list(bbox = matrix( rep(NA, 2), 2, 2, dimnames = list(NULL, c("min","max"))),
                                    proj4string = CRS(as.character(NA)),
                                    rings = list()),
    
    validity = function(object) {
        if (any(unlist(lapply(object@rings, function(x) !is(x, "Ring"))))) 
            stop("rings not Ring objects")
        if (any(duplicated(sapply(slot(object, "rings"), function(i) slot(i, "ID")))))
            return("non-unique Rings ID slot values")
        return(TRUE)
    }
)

#########################
###
### SpatialRingsDataFrame
###
#########################

setClass("SpatialRingsDataFrame",
    representation("SpatialRings", data = "data.frame"),
    validity = function(object) {
        if (!inherits(object@data, "data.frame"))
            stop("data should be of class data.frame")
        if (nrow(object@data) != length(object@rings))
            stop("number of rows in data.frame and SpatialRings don't match")
        return(TRUE)
    }
)

#########################
###
### SpatialCollections
###
#########################

setClassUnion("SpatialPointsNULL", c("SpatialPoints", "NULL")) 
setClassUnion("SpatialLinesNULL", c("SpatialLines", "NULL")) 
setClassUnion("SpatialRingsNULL", c("SpatialRings", "NULL")) 
setClassUnion("SpatialPolygonsNULL", c("SpatialPolygons", "NULL")) 

setClass("SpatialCollections",
    representation("Spatial", pointobj = "SpatialPointsNULL", 
							  lineobj = "SpatialLinesNULL",
							  ringobj = "SpatialRingsNULL",
							  polyobj = "SpatialPolygonsNULL",
							  plotOrder = "integer"),
    
    prototype = list(bbox = matrix( rep(NA, 2), 2, 2, dimnames = list(c("x","y"), c("min","max"))),
                     proj4string = CRS(as.character(NA)),
                     pointobj = NULL,
					 lineobj = NULL,
					 ringobj = NULL,
					 polyobj = NULL ),
    
    validity = function(object) {
		if (!is.null(object@pointobj)) validObject(object@pointobj)
		if (!is.null(object@lineobj)) validObject(object@lineobj)
		if (!is.null(object@ringobj)) validObject(object@ringobj)
		if (!is.null(object@polyobj)) validObject(object@polyobj)
		
		return(TRUE)
    }
)
