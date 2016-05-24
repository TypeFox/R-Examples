
setMethod("show","Cartesian2DCoordinateTransformation",
	function(object){
		cat("\n*** Object of the class '", class(object)[1], "' *** \n")
		nrowShow <- min(10,nrow(object@controlPoints))
		if(length(object@controlPoints) != 0){
			cat("\n* Slot 'controlPoints' (display limited to 10 rows): \n")
			colnames(object@controlPoints) = c("X source", "Y source", "X target", "Y target")
			print(object@controlPoints[1:nrowShow,],quote=FALSE)
		}else{
			cat("\n* Slot 'controlPoints' is empty \n")
		}
		if(length(object@parameters) == 0){
			cat("\n* Slot 'parameters' is empty \n");
		}else{
			names(object@parameters) = letters[seq( from = 1, to = length(object@parameters) )]
			cat("\n* Slot 'parameters': \n"); print (object@parameters)
		}
		cat("\n")
	}
)

### Getter for slot parameters
setGeneric("getParameters",function(object){standardGeneric ("getParameters")})
setMethod("getParameters","Cartesian2DCoordinateTransformation",
	function(object){
		if (length(object@parameters) == 0)
			return("Transformation parameters are unknown. You may need to call 'calculateParameters' first.")
		names(object@parameters) = letters[seq( from = 1, to = length(object@parameters) )]
		return(object@parameters)
	}
)

### Getter for slot residuals
setGeneric("getResiduals",function(object){standardGeneric ("getResiduals")})
setMethod("getResiduals","Cartesian2DCoordinateTransformation",
	function(object){
        if (nrow(object@residuals) == 0 && ncol(object@controlPoints) == 0)
   			return("Residuals cannot be calculated. You must provide control points first.")
		if (nrow(object@residuals) == 0 && ncol(object@controlPoints) > 0)
			return("Residuals are unknown. You may need to call 'calculateParameters' first.")
        if (object@residuals == 0 && length(object@parameters) > 0)
   			return("Residuals are all zero. There are not redundant control points to apply Least Squares.")
		return(object@residuals)
	}
)

### Getter for slot rmse
setGeneric("getRMSE",function(object){standardGeneric ("getRMSE")})
setMethod("getRMSE","Cartesian2DCoordinateTransformation",
	function(object){
        if (is.null(object@rmse) && ncol(object@controlPoints) == 0)
   			return("Root Mean Square Error cannot be calculated. You must provide control points first.")
		if (is.null(object@rmse) && ncol(object@controlPoints) > 0) 
			return("Root Mean Square Error is unknown. You may need to call 'calculateParameters' first.")
        if (object@rmse == 0 && length(object@parameters) > 0) 
   			return("Root Mean Square Error is zero. There are not redundant control points to apply Least Squares.")
		return(object@rmse)
	}
)

### Transform coordinates using calculated parameters
### Arguments:
### - coords is a matrix containing x and y values
### - object is either an "AffineTransformation" or "SimilarityTransformation" object
### Returns:
### A matrix containing transformed x and y values
###
"transformCoordinates" <- function(coords, object) {
    X <- object@parameters
    if (is(object, "AffineTransformation"))
        newCoords <- vapply(1:nrow(coords),
            FUN=function(x) c(c(coords[x,],1)%*%X[1:3],c(coords[x,],1)%*%X[4:6]),
            FUN.VALUE=c(0,0)
        )
    if (is(object, "SimilarityTransformation"))
	    newCoords <- vapply(1:nrow(coords),
            FUN=function(x) c(c(coords[x,],1)%*%X[1:3],c(coords[x,],1)%*%c(-X[2],X[1],X[4])),
            FUN.VALUE=c(0,0)
        )
	t(newCoords)
}

### Plot a transformed grid by using the transformation parameters
### Arguments:
### - object is either a "AffineTransformation" or "SimilarityTransformation" object
### - bbox is an SP bbox object, i.e. a 2x2 matrix with coordinates
### - numberOfPoints is the number of points to represent the grid
###
setGeneric("plotGridTransformation",function(object, bbox, numberOfPoints){
    standardGeneric ("plotGridTransformation")})
setMethod(f="plotGridTransformation",signature(object="Cartesian2DCoordinateTransformation"),
    definition=function(object, bbox, numberOfPoints){

        if (missing(object))
            stop("Please provide a transformation object as first argument.")
        if (missing(bbox))
            stop("Please provide a bounding box (bbox).")
        if (missing(numberOfPoints))
            stop("Please provide a number of points for representing the grid.")

        # Adapted from SP package (Class-Spatial.R)
        if (!is.matrix(bbox))
            stop("bbox should be a matrix")
        if (any(is.na(bbox)))
            stop("bbox should never contain NA values")
        if (any(!is.finite(bbox)))
            stop("bbox should never contain infinite values")
        if (any(bbox[,"max"] < bbox[,"min"]))
            stop("invalid bbox: max < min")
        # end of "Adapted from SP package"

        if (!isTRUE(all.equal(numberOfPoints,as.integer(numberOfPoints))))
            stop("numberOfPoints should be an integer")
        if (numberOfPoints <= 0)
            stop("numberOfPoints should be greater than zero")

		if (length(object@parameters) == 0)
			stop("Parameters have to be calculated before. Call 'calculateParameters' and try again.")

        offset=as.integer(sqrt(numberOfPoints))

        # Adapted from http://casoilresource.lawr.ucdavis.edu/drupal/node/433
        x <- seq(bbox[1,1],bbox[1,2], length.out=offset)
        y <- seq(bbox[2,1],bbox[2,2], length.out=offset)
        g = expand.grid(X=x, Y=y)
        ng = transformCoordinates(as.matrix(g), object)
        plot(g, cex=0.3, main='Transformed grid', col='red')
        points(ng, cex=.3, col='green')
        legend(bbox[1,1],bbox[2,2],legend=c("Original point", "Transformed point"),
            pch=c(16,16),pt.cex=.6,col=c('red', 'green'))
        # End of "Adapted from http://casoilresource.lawr.ucdavis.edu/drupal/node/433"
	}
)

### Apply the transformation to an SP object
### Arguments:
### - object is either a "AffineTransformation" or "SimilarityTransformation" object
### - sp.object is an object of type: SpatialPoints, SpatialPointsDataFrame, 
###      SpatialLines, SpatialLinesDataFrame, SpatialPolygons or SpatialPolygonsDataFrame
### Returns:
###   A transformed sp.object
###
setGeneric("applyTransformation",function(object, sp.object){
    standardGeneric ("applyTransformation")})
setMethod(f="applyTransformation",signature(object="Cartesian2DCoordinateTransformation"),
    definition=function(object, sp.object){

        if (missing(object))
            stop("Please provide a transformation object as first argument.")
        if (missing(sp.object))
            stop("Please provide an SP object as second argument.")

        if (!class(sp.object) %in% c('SpatialPoints', 'SpatialPointsDataFrame', 
            'SpatialLines', 'SpatialLinesDataFrame', 'SpatialPolygons', 
            'SpatialPolygonsDataFrame'))
            stop('Transformation can be applied on objects of type SpatialPoints, SpatialPointsDataFrame, SpatialLines, SpatialLinesDataFrame, SpatialPolygons or SpatialPolygonsDataFrame.')

		if (length(object@parameters) == 0)
			stop("Parameters have to be calculated before. Call 'calculateParameters' and try again.")

        if (!is.na(proj4string(sp.object)) && !is.projected(sp.object))
            stop("The SP object cannot have a geographic Coordinate Reference System (CRS).")


        bDataFrame = ifelse(class(sp.object) %in% c('SpatialPointsDataFrame', 'SpatialLinesDataFrame', 
            'SpatialPolygonsDataFrame'), TRUE, FALSE)
            
        if (bDataFrame)
            df = sp.object@data

        rs = CRS(proj4string(sp.object))

        # Do transform!
        if (is(sp.object, 'SpatialPoints')){
            newCoords = transformCoordinates(coordinates(sp.object), object)
    
            if (bDataFrame){
                newSPObject=SpatialPointsDataFrame(coords=newCoords, data=df, proj4string=rs)
            }else{
                newSPObject=SpatialPoints(coords=newCoords, proj4string=rs)
            }
        }

        if (is(sp.object, 'SpatialLines')){
            newLines = lapply( sp.object@lines, 
            	function(objLines) Lines(lapply( objLines@Lines,  
            		function(objLine) Line(transformCoordinates(objLine@coords, object))),
            		ID=objLines@ID)
            )
            newSPObject = SpatialLines(newLines, proj4string=rs)

            if (bDataFrame)
                newSPObject = SpatialLinesDataFrame(newSPObject, data=df)
        }

        if (is(sp.object, 'SpatialPolygons')){
            newPolygons = lapply( sp.object@polygons, 
            	function(objPolygons) Polygons(lapply( objPolygons@Polygons,  
            		function(objPolygon) Polygon(transformCoordinates(objPolygon@coords, object), hole=objPolygon@hole)),
            		ID=objPolygons@ID)
            )
            newSPObject = SpatialPolygons(newPolygons, pO=sp.object@plotOrder, proj4string=rs)

            if (bDataFrame)
                newSPObject = SpatialPolygonsDataFrame(newSPObject, data=df)
        }

		return(newSPObject)
	}
)

