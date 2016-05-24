setClass("SpatialGDAL",
    representation("Spatial", grid = "GridTopology", grod = "GDALReadOnlyDataset", name = "character"),
    validity = function(object) {
        if (is.null(object@grod))
            stop("grod is NULL; this should not happen")
        return(TRUE)
    }
)
setClass("SpatialGDALWrite", "SpatialGDAL")

open.SpatialGDAL = function(con, ..., silent = FALSE) {
	if (nchar(con) == 0) stop("empty file name")

	grod = GDAL.open(con, read.only = TRUE, silent=silent)

	d = dim(grod)
	if (!silent) {
		cat(paste(con, "has GDAL driver", getDriverName(getDriver(grod)),"\n"))
		cat(paste("and has", d[1], "rows and", d[2], "columns\n"))
	}
#	p4s <- .Call("RGDAL_GetProjectionRef", grod, PACKAGE="rgdal")
	p4s <- getProjectionRef(grod, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=NULL)
	if (nchar(p4s) == 0) p4s <- as.character(NA)
	gt = .Call('RGDAL_GetGeoTransform', grod, PACKAGE="rgdal")
        if (attr(gt, "CE_Failure")) warning("GeoTransform values not available")
	# [1] 178400     40      0 334000      0    -40
	if (any(gt[c(3,5)] != 0.0)) 
		stop("rotated grid cannot be read; try readGDAL to read as points")
	cellsize = abs(c(gt[2],gt[6]))
	ysign <- sign(gt[6])
	output.dim <- dim(grod)[1:2] # rows=nx, cols=ny
	half.cell <- c(0.5, 0.5)
	co.x <- gt[1] + (half.cell[2]) * cellsize[1]
	co.y <- ifelse(ysign < 0, gt[4] + (ysign*((output.dim[1]) + (ysign*half.cell[1]))) * abs(cellsize[2]),
		gt[4] + (ysign*((offset[1]) + (ysign*half.cell[1]))) * 
		abs(cellsize[2]))
	cellcentre.offset <- c(x=co.x, y=co.y)
	grid = GridTopology(cellcentre.offset, cellsize, rev(output.dim))
	#bbox = bbox(SpatialGrid(grid)),

	x = SpatialPoints(rbind(grid@cellcentre.offset, 
			grid@cellcentre.offset + (grid@cells.dim - 1) * grid@cellsize))
	x@bbox[,1] = x@bbox[,1] - 0.5 * grid@cellsize
	x@bbox[,2] = x@bbox[,2] + 0.5 * grid@cellsize

	data = new("SpatialGDAL", 
		bbox = x@bbox,
		proj4string = CRS(p4s), 
		grid = grid,
		# data = data.frame(), 
		grod = grod,
		name = con)
	return(data)
}

close.SpatialGDAL = function(con, ...) {
	GDAL.close(con@grod)
	invisible(NULL)
}
close.SpatialGDALWrite = close.SpatialGDAL 

copy.SpatialGDAL = function(dataset, fname, driver = 
		getDriver(dataset@grod), strict = FALSE, options = NULL, silent = FALSE)
{
	if (nchar(fname) == 0) 
		stop("empty file name")
	assertClass(dataset, 'SpatialGDAL')
	# grod = GDAL.open(con, read.only = FALSE)

	if (is.character(driver)) 
		driver <- new("GDALDriver", driver)
	if (nchar(fname) == 0) 
		stop("empty file name")

	if (!is.null(options) && !is.character(options))
		stop("options not character")

	grod <- new('GDALTransientDataset',
		handle = .Call('RGDAL_CopyDataset',
			dataset@grod, driver, as.integer(strict),
			as.character(options), fname, PACKAGE="rgdal"))

	data = new("SpatialGDALWrite", 
		bbox = bbox(dataset),
		proj4string = CRS(proj4string(dataset)), 
		grid = dataset@grid, grod = grod, name = fname)
	return(data)
}

setAs("SpatialGDAL", "SpatialGridDataFrame",
	function(from) { 
#		print("doing the coerce...")
		from[]
	}
)

setAs("SpatialGDAL", "SpatialPixelsDataFrame",
	function(from) as(from[], "SpatialPixelsDataFrame")
)

setMethod("[", "SpatialGDAL",
	function(x, i, j, ... , drop = FALSE)
		x@grod[i = i, j = j, ...]
)

setMethod("[[", c("SpatialGDAL", "ANY", "missing"),
    function(x, i, j, ...)
        x[,,i][[1]] # efficient: reads only band i
)

# avoid inheritance:
setMethod("$", "SpatialGDAL",
    function(x, name)
        stop("use [[ with numeric index to select single band")
)

# avoid inheritance:
setReplaceMethod("$", "SpatialGDAL",
    function(x, name, value)
		stop("no replacement method available")
)

# avoid inheritance:
setReplaceMethod("[[", c("SpatialGDAL", "ANY", "missing", "ANY"),
    function(x, i, j, value)
		stop("no replacement method available")
)

setMethod("summary", "SpatialGDAL",
	function(object, ...) {
		obj = list()
		obj$class = class(object)
		obj$grid = object@grid
		obj$name = object@name
		class(obj) = "summary.SpatialGDAL"
		obj
	}
)

print.summary.SpatialGDAL = function(x, ...) {
	cat(paste("object of class", x$class, "\n"))
	cat(paste("file name:", x$name, "\n"))
	print(x$grid)
}

setReplaceMethod("[", "SpatialGDALWrite", function(x, i, j, ..., value) {
	ncol = gridparameters(x@grid)[1,"cells.dim"]
	nrow = gridparameters(x@grid)[2,"cells.dim"]
	if (!is.numeric(value)) stop("Numeric bands required")
	if (missing(i)) i = 1:nrow
	if (missing(j)) j = 1:ncol
	i = i - 1 # y/row index -- zero offset
	j = j - 1 # x/cols -- zero offset
	if (length(i) * length(j) != length(value))
		stop("lengths do not match")
	dots = list(...)
	if (length(dots) > 0L)
		band = dots[[1]]
	else
		band = 1
#	mvFlag = .Call("RGDAL_GetNoDataValue", x@grod)
#	if (!is.na(mvFlag))
#		data[is.na(data)] = mvFlag
	if (!is.matrix(value))
		value = matrix(value, length(i), length(j))
	putRasterData(x@grod, value, band, c(i[1], j[1])) # offset y, x
	x
})
