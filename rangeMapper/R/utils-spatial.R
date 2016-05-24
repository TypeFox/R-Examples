#' Convert WKT polygons to SpatialPolygonsDataFrame
#'
#' Convert a data.frame containing WKT polygons to a \code{SpatialPolygonsDataFrame}.
#'
#'
#' @param dat   \code{data.frame}
#' @param geom  is the name (character vector) of the column in the
#'              \code{data.frame} containing the geometry.
#' @param id    is the name (character vector) of the column in the
#'              \code{data.frame} identifying the polygon.  when \code{id} is not unique
#'              then polygons are combined using \code{\link[rgeos]{gUnionCascaded}}.
#' @return      a \code{\link[sp]{SpatialPolygonsDataFrame} }object.
#' @export
#' @examples
#'
#' require(rangeMapper)
#' require(rgeos)
#'
#' # generate a few random polygons
#' randPoly = function(mean, sd) {
#'   writeWKT(
#'     gConvexHull(
#'      readWKT(paste("MULTIPOINT (",
#'              paste(apply(matrix(rnorm(n= 100, mean, sd), ncol = 2), 1,
#'              paste, collapse = ' '), collapse = ","), ")"))))
#' }
#' n = 50
#' d = data.frame( nam = sample(letters, n, TRUE),
#'                range = mapply(randPoly, mean = sample(1:2, n, TRUE),
#'                sd = sample(1:2/5, n, TRUE) ))
#'
#'
#' X = WKT2SpatialPolygonsDataFrame(d, 'range', 'nam')
#'
#'
#' dbcon = rangeMap.start(file = "test.sqlite", overwrite = TRUE, dir = tempdir() )
#' global.bbox.save(con = dbcon, bbox = X)
#' gridSize.save(dbcon)
#' canvas.save(dbcon)
#' processRanges(spdf = X, con =  dbcon, ID = "nam")
#' rangeMap.save(dbcon)
#' plot(rangeMap.fetch(dbcon))
#'
WKT2SpatialPolygonsDataFrame <- function(dat, geom, id) {
	dl = split(dat, dat[, id])

	o = lapply(dl, function(x) {

		p = mapply(readWKT, text = x[, geom], id = 1:nrow(x), USE.NAMES = FALSE )
			if(length(p) == 1) {
				p = p[[1]]
				p = spChFIDs(p, as.character(x[1, id]))
				}

			if(length(p) > 1) {
				p = do.call(sp::rbind.SpatialPolygons, p)
				p = gUnionCascaded(p, id = as.character(x[1, id]) )
				}
	p
	})

	X = do.call(sp::rbind.SpatialPolygons, o)
	dat = data.frame(id = sapply(slot(X, "polygons"), function(x) slot(x, "ID")) )
	row.names(dat ) = dat$id
	names(dat) = id
	X = SpatialPolygonsDataFrame(X, data =  dat)
	X

	}

#' Vertices of a SpatialPolygonsDataFrame
#'
#' Extract vertices from a \link[sp]{SpatialPolygonsDataFrame} and optionally
#' applies an aggregating function to each Polygon.
#'
#' @param 	 object  An object.
#' @param 	 FUN  A function.
#' @return   A \link[sp]{SpatialPointsDataFrame} containing an id column
#'           corresponding to each extracted Polygon.
#' @export
#' @examples
#' require(rangeMapper)
#' require(rgdal)
#' f = system.file(package = "rangeMapper", "extdata", "wrens", "vector")
#' # path to Campylorhynchus_gularis breeding range:
#' camgul = selectShpFiles(f, ogr = TRUE, polygons.only = TRUE)[6, ]
#' r = readOGR(camgul$dsn, camgul$layer)
#' mp = vertices(r, mean)
#' v = vertices(r)
#'
#' plot(r)
#' points(mp, col = 2, pch = 3, cex = 2)
#' points(v, pch = 3, cex = .5)
#'
#' @rdname WKT2SpatialPolygonsDataFrame
setGeneric("vertices", function(object, FUN)  standardGeneric("vertices") )

#' @rdname WKT2SpatialPolygonsDataFrame
setMethod("vertices", "SpatialPolygons",
	function(object, FUN) {
		d = lapply( unlist(lapply(slot(object, "polygons"), function(P) slot(P, "Polygons"))), function(cr) slot(cr, "coords") )
		d = lapply(d, function(x) { dimnames(x)[[2]] = c('x', 'y'); x} )
		d = lapply(d, function(x) x[-nrow(x), , drop = FALSE])
		d = mapply("cbind", 1:length(d), d, SIMPLIFY = FALSE)
		if(!missing(FUN))
			d = lapply(d, function(x) apply(x ,2, FUN) )

		d = data.frame(do.call("rbind", d))

		coordinates(d) = ~ x+y
		proj4string(d) = CRS(proj4string(object))
		names(d) = "id"
		d
	})

#' A container of functions to apply on a \code{SpatialPolygons} object
#'
#' This is a convenience function returning a named \code{list} of
#' functions.
#'
#' The function returns a named list so any additional functions should be
#' given as rangeTraits(funName1 = FUN1, funName2 = FUN2) where FUN1, FUN2 are
#' \code{\link{SpatialPolygons}} extractor functions.
#'
#' @param \dots       functions, given as myfun = FUN, to apply on a
#'                    \code{\link{SpatialPolygons}} object
#' @param use.default If \code{TRUE}, the default, the output list contains
#'                    functions to extract Area, Median, Min and Max extent of the
#'                    \code{\link{SpatialPolygons}} object. This option is ignored
#'                    if no functions are given.
#' @return            Returns a named list containing extractor functions to apply on
#'                    \code{\link[sp]{SpatialPolygons}} objects.
#' @seealso           \code{\link{processRanges}}.
#' @export
#' @examples
#' summary(rangeTraits(use.default = FALSE))
#'
#' f = system.file(package = "rangeMapper", "extdata", "wrens", "vector")
#' troaed = selectShpFiles(f, ogr = TRUE,
#' 	polygons.only = TRUE)[71, ] # path to Troglodytes_aedon
#' require(rgdal)
#' r = readOGR(troaed$dsn, troaed$layer)
#'
#' # Beware of the value returned for Area!
#' sapply(rangeTraits(), function(x) x(r) )
#'
#' # Define an extra function to compute correct Area
#' Area2 = function(x) {
#' x = spTransform(x,
#' CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
#' 	)
#'
#' sum(sapply(slot(x, "polygons"), function(x) slot(x, "area") ))
#' }
#'
#' sapply(rangeTraits(Area_sqm = Area2), function(x) x(r) )
#'
#'
#'
rangeTraits <- function(..., use.default = TRUE) {

	Area     = function(spdf) sum(sapply(slot(spdf, "polygons"), function(x) slot(x, "area") ))
	Median_x = function(spdf) median(coordinates(vertices(spdf))[, 1])
	Median_y = function(spdf) median(coordinates(vertices(spdf))[, 2])
	Min_x    = function(spdf) min(coordinates(vertices(spdf))[, 1])
	Min_y    = function(spdf) min(coordinates(vertices(spdf))[, 2])
	Max_x    = function(spdf) max(coordinates(vertices(spdf))[, 1])
	Max_y    = function(spdf) max(coordinates(vertices(spdf))[, 2])


	res = list(Area = Area, Median_x = Median_x, Median_y = Median_y, Min_x = Min_x, Min_y = Min_y, Max_x = Max_x, Max_y = Max_y)

	x = list(...)
	if(length(x) > 0) {
		 if(length(names(x)) != length(x)) stop (dQuote("..."), " elements should be named, e.g. myFun = abc")
		 if( !all(sapply(x, is.function))) stop (dQuote("..."), " elements should be functions.")
		 if(use.default) res = c(res, x)
	}
	res
	}

#' rangeOverlay
#'
#' @param spp     a SpatialPolygons* object.
#' @param canvas  a SpatialPointsDataFrame.
#' @param name    a character vector (a bioid)
#'
#' @return        a data.frame with two columns: id (canvas id) and bioid.
#' @export
#'
rangeOverlay <- function(spp, canvas, name) {

	if(inherits(spp, 'SpatialPolygonsDataFrame'))
		spp = as(spp, 'SpatialPolygons')

	overlayRes = which(!is.na(over(canvas, spp)))


	if(length(overlayRes) > 0) { 	# do grid interpolation
		sp = canvas[overlayRes, ]
		o = data.frame(id = sp$id, bioid = rep(name, nrow(sp)) )
		}

	if(length(overlayRes) == 0) { 	# the polygons are smaller than the grid cells:  snap to the nearest points
			xy = vertices(spp, FUN = mean)
			nn = spDists(canvas, xy)
			mins = apply(nn, 2, min)
			res = vector(mode = 'numeric', length = length(mins))
			for(i in 1:length(res)) {
				res[i] = which(nn[,i] == mins[i])
				}
			res = unique(res)

			sp = canvas[res, ]

			o = data.frame(id = sp$id, bioid = rep(name, nrow(sp@coords)) )
		}

	return(o)
	}


# undocumented functions
extract.p4s <- function(ShpFiles) {
	#ShpFiles = selectShpFiles(paste(system.file(package="rangeMapper"), "extdata", "wrens", "vector", sep = .Platform$file.sep))
		fl = split(ShpFiles, ShpFiles$layer)
		unlist(lapply(fl, FUN = function(x) OGRSpatialRef(x[,1], x[,2])  ))
	}

rect2spp <- function(xmin, xmax, ymin, ymax) {
	bb = cbind(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin) )
	SpatialPolygons(Srl = list(Polygons(list(Polygon(bb)), "bb")) )
	}

proj4string_is_identical <- function(a, b){
	identical(gsub(" ", "", a), gsub(" ", "", b) )
	}

rmap.frame_to_SpatialPixelsRangeMap <- function(map, proj4string, names){
	map = new("SpatialPixelsRangeMap",
			SpatialPixelsDataFrame(
				points = map[, c('x', 'y'), with = FALSE] %>% as.matrix,
				data   = map[,  setdiff(names(map), c('x', 'y')), with = FALSE]  %>% as.data.frame,
				proj4string = CRS(proj4string) ),
			mapvar = names)
	}