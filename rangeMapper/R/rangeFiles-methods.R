setGeneric("rangeFiles", function(object, ...)   					standardGeneric("rangeFiles") )

setMethod("rangeFiles",
	signature = "rangeFiles",
	definition = function(object, ...){
	dir.list = list.files(object@dir, recursive = TRUE, full.names = TRUE, pattern = ".shp$")

	if(object@polygons.only){
	nfo = lapply(dir.list, maptools::getinfo.shape)
		is.poly =  unlist(lapply(nfo, function(x) x$type == 5))
		dir.list = dir.list[is.poly]
		}

	if(object@ogr)
		res = data.frame(dsn = dirname(dir.list), layer = gsub(".shp", "", basename(dir.list)), stringsAsFactors = FALSE) else
		res = dir.list
		res
		} )

#' Select (recursively) shape files
#'
#' Returns the file path to all \sQuote{.shp} polygons in a directory.
#'
#' @param dir    character string specifying the directory containing .shp files.
#' @param \dots  currently ignored
#' @return       Either a \code{\link{data.frame}} or a character vector is returned.
#' @note         The function uses \code{\link[maptools]{getinfo.shape}} to only select
#'               polygon files (aka type 5).
#' @examples
#' f = system.file(package="rangeMapper", "extdata", "wrens", "vector")
#' res = selectShpFiles(f, ogr = TRUE, polygons.only = TRUE)
#' head(res)
#'
#' @export selectShpFiles
selectShpFiles <- function(dir,  ...) {
	rangeFiles(new("rangeFiles", dir =dir, ...))
	}


