readShapeSpatial <- function(fn, proj4string=CRS(as.character(NA)), 
	verbose=FALSE, repair=FALSE, IDvar=NULL, force_ring=FALSE, 
	delete_null_obj=FALSE, retrieve_ABS_null=FALSE) {
	shinfo <- getinfo.shape(fn)
	if (verbose) print(shinfo)
	type <- shinfo[[2]]
	types <- c("Point", NA, "PolyLine", NA, "Polygon", NA, NA, 
	    "MultiPoint", NA, NA, "PointZ", NA, "PolyLineZ", NA, 
	    "PolygonZ", NA, NA, "MultiPointZ", NA, NA, "PointM", NA, 
	    "PolyLineM", NA, "PolygonM", NA, NA, "MultiPointM", NA, NA, 
	    "MultiPatch")
	typeSh <- types[type]
	if (typeSh == "Point" || typeSh == "PointZ" || typeSh == "MultiPoint") {
	    res <- readShapePoints(fn=fn, proj4string=proj4string, 
		verbose=verbose, repair=repair)
	} else if (typeSh == "PolyLine" || typeSh == "PolyLineZ") {
	    res <- readShapeLines(fn=fn, proj4string=proj4string, 
		verbose=verbose, repair=repair)
	} else if (typeSh == "Polygon" || typeSh == "PolygonZ") {
	    res <- readShapePoly(fn=fn, IDvar=IDvar, proj4string=proj4string, 
		verbose=verbose, repair=repair, force_ring=force_ring, 
		delete_null_obj=delete_null_obj, 
		retrieve_ABS_null=retrieve_ABS_null)	    
	} else stop("File type cannot be read")
	res
}

writeSpatialShape  <- function(x, fn, factor2char = TRUE, max_nchar=254) {
	if (is(x, "SpatialPolygonsDataFrame")) {
	    writePolyShape(x=x, fn=fn, factor2char=factor2char, 
		max_nchar=max_nchar)
	} else if (is(x, "SpatialLinesDataFrame")) {
	    writeLinesShape(x=x, fn=fn, factor2char=factor2char, 
		max_nchar=max_nchar)
	} else if (is(x, "SpatialPointsDataFrame")) {
	    writePointsShape(x=x, fn=fn, factor2char=factor2char, 
		max_nchar=max_nchar)
	} else {
            stop("x is a", class(x), "object, not a compatible Spatial*DataFrame")
        }
}
