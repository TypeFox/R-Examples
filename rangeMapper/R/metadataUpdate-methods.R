setGeneric("metadataUpdate", function(object, FUN,name, map, parallel, ...)  standardGeneric("metadataUpdate") )

# Method 1: using  a SpatialGridDataFrame object
setMethod("metadataUpdate",
		signature = c(object = 'rangeMap', FUN = 'function', name = 'character',
					 map = 'SpatialGridDataFrame', parallel = 'missing'),
		definition = function(object, FUN, name, map, ...){ # ... goes to FUN
			Startprocess = Sys.time()

			# only use 1st band of the raster
			if(length(names(map)) > 1) {
				warning( paste("The SpatialGridDataFrame has more than one band, only", dQuote(names(map)[1]), "will be used!") )
			map = map[names(map)[1]]
			}

			# check if metadata_ranges is populated
			empty = dbGetQuery(object@CON, "SELECT count(bioid) from metadata_ranges") == 0
			if(empty)
				dbGetQuery(object@CON, "INSERT INTO metadata_ranges SELECT distinct bioid from ranges")

			# add new column
			if( is.element(name, c(names(dbGetQuery(object@CON, 'SELECT * FROM metadata_ranges limit 1')), SQLKeywords(object@CON) ) ) )
				stop(paste(dQuote(name), "allready exists or is not a valid SQL name!") )
			name = RSQLite::make.db.names(object@CON, name)
			dbGetQuery(object@CON, paste("ALTER TABLE metadata_ranges ADD COLUMN", name ,"NUMERIC;") )

			# get ranges listing
			mr = dbGetQuery(object@CON, "SELECT bioid from metadata_ranges")$bioid

			# loop through ranges,  apply FUN and update metadata_ranges

			message("Updating metadata_ranges ...")

			for( i in 1:length(mr) ) {
				idi = mr[i]
				ri = rangeFetch(object, idi)
				sraster = map[!is.na(over(map, ri)), ] #  SpatialGridDataFrame subset
				x = as.numeric(sraster@data[,1])
				res = FUN(x, ...) #apply FUN

				if( any(length(res) > 1 | res%in%c(-Inf, Inf)) ) {
					dbRemoveField(object@CON, 'metadata_ranges', name)
					stop( paste("FUN returned ", dQuote(res), ". It should only return a finite numeric vector of length 1.", sep = '') )
				}

				if(!is.na(res))
					dbGetQuery(object@CON, paste('UPDATE metadata_ranges SET' ,name, '=' ,res, 'WHERE bioid =' ,shQuote(idi) ))

			}

		# last Msg
		message( paste("Done in ", round(difftime(Sys.time(), Startprocess, units = "mins"),1), "mins") )
		} )

#' Updates metadata table
#'
#' Updates \code{metadata_table} of a \code{rangeMapper} project \emph{after}
#' importing ranges with \code{\link{processRanges}}.
#'
#' @param rangeMap  A \code{\link{rangeMap}} object.
#' @param FUN       Function used to aggregate the map values corresponding to each
#'                  range
#' @param name      The name of the new \code{metadata_table} field containing the
#'                  variable computed by \code{FUN}
#' @param map       Single-band \code{\link[sp]{SpatialGridDataFrame}} object
#' @param overwrite If set to \code{TRUE} the the values of the field are
#' replaced
#' @param \dots     extra arguments (e.g. \code{na.rm = TRUE}) to be passed to FUN.
#' @return          NULL.
#' @note            In order to compute taxa-level metadata which are not dependent on the
#'                  project's resolution use \code{\link{processRanges}} with a \code{metadata}
#'                  argument. See \code{\link{rangeTraits}} for more details. \cr The method can
#'                  be extended to work with raster or vector objects (e.g. lines, polygons,
#'                  points) using overlaying functions in the package \code{raster} and
#'                  \code{rgeos} respectively.
#' @export
#' @examples
#'
#' require(rangeMapper)
#' require(rgdal)
#' # data
#' spdf = readOGR(system.file(package = "rangeMapper",
#' 	"extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)
#' rloc = system.file(package = "rangeMapper", "extdata",
#' 	"etopo1", "etopo1_Americas.tif")
#' r = readGDAL(rloc, output.dim = c(50, 50))
#' spdf = spTransform(spdf, CRS(proj4string(r)) )
#'
#' # the project
#' dbcon = rangeMap.start(file = "wrens.sqlite", overwrite = TRUE,
#' 	dir = tempdir() )
#' rmap = new("rangeMap", CON = dbcon)
#' global.bbox.save(con = dbcon, bbox = spdf )
#' gridSize.save(dbcon, gridSize = 300000)
#' canvas.save(dbcon)
#' processRanges(spdf = spdf, con =  dbcon, ID = "sci_name" )
#'
#' # metadata.update
#' metadata.update (rmap,
#'  FUN = function(x, ...) {
#'  	res = diff(range(x, ...))
#'  	if( !is.finite(res)) res = 0
#'  	res
#'  	},
#' 	name = 'AltitudeRange', map = r, na.rm = TRUE, overwrite = TRUE)
#' # plot
#' mr = RSQLite::dbGetQuery(dbcon, 'select * from metadata_ranges')
#' maxRangeSp = mr[mr$AltitudeRange== max(mr$AltitudeRange), 'bioid']
#' image(r)
#' plot(rangeFetch(rmap, maxRangeSp), add = TRUE, border = 4, lwd = 3)
#' title(main = maxRangeSp)
#'
metadata.update  <- function(rangeMap, FUN, name, map, overwrite = FALSE,...){
	if(overwrite)
	dbRemoveField(rangeMap@CON, 'metadata_ranges', name)
	metadataUpdate(object = rangeMap, name = name, FUN = FUN, map = map, ...)
   }








