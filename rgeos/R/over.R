order_relations = function(rel, minDimension) {
	stopifnot(minDimension %in% 0:2)
	rel = sapply(rel, function(x)
			paste0(substring(x, c(1,4),c(2,5)), collapse=""))
		# our interest is in chars
		# 1-2 = inner of x with inner/border of y
		# 4-5 = border of x with inner/border of y
	ret = vector("numeric", length(rel)) * NA
	for (d in minDimension:2) {
		r = regexpr(as.character(d), rel, fixed = TRUE)
		sel = which(r != -1)
		if (length(sel) > 0)
			ret[sel] = 4 - r[sel] + 4 * d
	}
	order(ret, decreasing = TRUE, na.last = NA)
}

listifyMatrix = function(x) { # put columns in list elements
	if (!is.list(x)) {
		if (length(x) == 0)
			return(list(x))
		if (!is.matrix(x)) { # vector!
			nm = names(x)
			x = matrix(x, 1, length(x))
		} else
			nm = dimnames(x)[[2]]
		x = lapply(1:ncol(x), function(i) x[,i])
		names(x) = nm
	}
	x
}

overGeomGeom = function(x, y, returnList = FALSE, fn = NULL, ..., minDimension = -1) {
	stopifnot(identicalCRS(x, y))
	if (gridded(x))
		x = as(x, "SpatialPolygons")
	if (gridded(y))
		y = as(y, "SpatialPolygons")

	if (minDimension %in% 0:2)
		ret = apply(gRelate(x, y, byid = TRUE), 2, order_relations, minDimension = minDimension)
	else
		ret = apply(gIntersects(x, y, byid = TRUE), 2, which)
	ret = listifyMatrix(ret) # if not already list, create one now
	if (! returnList) # pick first, or NA if length is 0:
		sapply(ret, function(x) (x)[1])
	else 
		ret
}

# taken from: overDFGeneric in sp; 
# if modified here, consider modifying there as well!
overGeomGeomDF = function(x, y, returnList = FALSE, fn = NULL, ..., minDimension = -1) {
    r = overGeomGeom(x, y, returnList = TRUE, minDimension = minDimension)
    #ret = sp:::.overDF(r, y@data, length(x), returnList, fn, ...)
	#  length(x) differs from length(r) in case of SpatialMultiPoints!!!
	#  reason to change is sp::overMultiPoints
    ret = overDF_for_rgeos(r, y@data, length(r), returnList, fn, ...)
    if (!returnList)
        row.names(ret) = row.names(r)
    ret
}

#setMethod("over",
#    signature(x = "SpatialPoints", y = "SpatialPolygons"),
#	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialPoints", y = "SpatialLines"),
	        overGeomGeom)
#setMethod("over",
#    signature(x = "SpatialPoints", y = "SpatialPoints"),
#	        overGeomGeom)
#setMethod("over",
#    signature(x = "SpatialPolygons", y = "SpatialPoints"),
#	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialPolygons", y = "SpatialLines"),
	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialPolygons", y = "SpatialPolygons"),
	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialPoints"),
	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialPolygons"),
	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialLines"),
	        overGeomGeom)

# all with DataFrame:
#setMethod("over",
#    signature(x = "SpatialPoints", y = "SpatialPolygonsDataFrame"),
#	        overGeomGeomDF)
setMethod("over",
    signature(x = "SpatialPoints", y = "SpatialLinesDataFrame"),
	        overGeomGeomDF)
#setMethod("over",
#    signature(x = "SpatialPoints", y = "SpatialPointsDataFrame"),
#	        overGeomGeomDF)
#setMethod("over",
#    signature(x = "SpatialPolygons", y = "SpatialPointsDataFrame"),
#	        overGeomGeomDF)
setMethod("over",
    signature(x = "SpatialPolygons", y = "SpatialLinesDataFrame"),
	        overGeomGeomDF)
setMethod("over",
    signature(x = "SpatialPolygons", y = "SpatialPolygonsDataFrame"),
	        overGeomGeomDF)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialPointsDataFrame"),
	        overGeomGeomDF)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialPolygonsDataFrame"),
	        overGeomGeomDF)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialLinesDataFrame"),
	        overGeomGeomDF)

# lines & grids:
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialPixels"),
	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialGrid"),
	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialPixelsDataFrame"),
	        overGeomGeomDF)
setMethod("over",
    signature(x = "SpatialLines", y = "SpatialGridDataFrame"),
	        overGeomGeomDF)
setMethod("over",
    signature(x = "SpatialPixels", y = "SpatialLines"),
	        overGeomGeom)
setMethod("over",
    signature(x = "SpatialGrid", y = "SpatialLinesDataFrame"),
	        overGeomGeomDF)
