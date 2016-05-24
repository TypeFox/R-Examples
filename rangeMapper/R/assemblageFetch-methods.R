
#' @include projectINI-methods.R
NULL

#' Retrieves the species set of an arbitrary canvas cell
#'
#' \code{assemblageFetch} retrieves the species set of an arbitrary canvas cell
#' optionally with the associated life history data
#'
#' @param object   A \code{connection} object.
#' @param xy       A \code{\link[sp]{SpatialPoints}} object.
#' @param BIO      The name of the \code{BIO_table} containing species life-history data.
#' @return         A \code{data.frame} containing the bioid (e.g. species names),
#'                 the canvas id and optionally any associated life history data contained
#'                 in the \code{BIO_table} table.
#' @export
#' @examples
#'
#' require(rangeMapper)
#' require(rgdal)
#'
#' projName = "wrens.sqlite"
#' projLoc = paste(tempdir(), projName, sep = .Platform$file.sep)
#'
#' dbcon = rangeMap.start(file = projName,dir = tempdir() , overwrite = TRUE)
#' f = system.file(package = "rangeMapper", "extdata", "wrens", "vector_combined")
#' r = readOGR(f, "wrens", verbose = FALSE)
#' global.bbox.save(con = dbcon, bbox = r)
#' gridSize.save(dbcon, gridSize = 3)
#' canvas.save(dbcon)
#' data(wrens)
#' bio.save(con = dbcon, loc = wrens ,  ID = "sci_name")
#' processRanges(spdf = r, con =  dbcon, ID = "sci_name")
#' rangeMap.save(dbcon)
#'
#' sr = rangeMap.fetch(dbcon)
#' image(sr, axes = TRUE); grid()
#'
#' p = list(x = -76.39, y = 9.26)
#' # or use locator:  p =  locator(1)
#'
#' xy = SpatialPoints( do.call(cbind, p), proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs ") )
#' af = assemblageFetch(rangeMap(projLoc) , xy)
#' points(p, col = 4, cex = 2)
#' print(af)
#'
#' af = assemblageFetch(rangeMap(projLoc) , xy, "wrens")
#' print(af[, c(1, 4, 6:8)])
#'
#'
setGeneric("assemblageFetch", function(object, xy, BIO)		standardGeneric("assemblageFetch") )

#' @rdname  assemblageFetch

setMethod("assemblageFetch",
	signature  = c(object = "rangeMap", xy = "SpatialPoints", BIO = "missing"),
	definition = function(object, xy) {

		# CANVAS
		cnv = canvas.fetch(object@CON)

		#Assembladge IDs
		assembl_id = over(xy, cnv)$id

		# buid sql
		sql = paste("SELECT * FROM ranges WHERE id in(",	paste(assembl_id, collapse = ",")  ,")")

		#fetch assambladges
		A = dbGetQuery(object@CON, sql)

		return(A)
		}
   )

#' @rdname  assemblageFetch

setMethod("assemblageFetch",
	signature  = c(object = "rangeMap", xy = "SpatialPoints", BIO = "character"),
	definition = function(object, xy, BIO) {

		# CANVAS
		cnv = canvas.fetch(object@CON)

		# BIO_table
		biotabs = dbGetQuery(object@CON, "SELECT * FROM sqlite_master WHERE type='table' and name like 'BIO_%' ")$name
		if(BIO%in%biotabs) stop(paste(dQuote(BIO), "is not a BIO_table"))
		BIO = paste("BIO", BIO, sep = "_")
		biotab_id = extract.indexed(object@CON, BIO)

		#Assembladge IDs
		assembl_id = over(xy, cnv)$id

		# buid sql
		sql = paste("SELECT * FROM ranges r LEFT JOIN", BIO, "ON", paste(BIO,biotab_id, sep = ".") , "= r.bioid",
						"WHERE r.id in(",	paste(assembl_id, collapse = ",")  ,")")

		#fetch assambladges
		A = dbGetQuery(object@CON, sql)
		A$bioid = NULL

		return(A)
		}
 )

