setGeneric("bioSave", function(object, ...)  	standardGeneric("bioSave") )

setMethod("bioSave",
	signature  = "bioSaveFile",
		definition = function(object) {

		tableName = paste(object@BIO, object@tableName, sep = "")

		d = read.table(object@loc, sep = object@sep, header = TRUE, stringsAsFactors = FALSE)

		nam = d[, object@ID]

		ranges.nam = dbGetQuery(object@CON, "select distinct bioid from ranges")$bioid

		d$has_range= is.element(nam, ranges.nam)

		res = dbWriteTable(object@CON ,tableName , d, row.names = FALSE)

		if(res) {
			dbGetQuery(object@CON,(paste("CREATE  INDEX", paste(tableName, object@ID, sep = "_") , "ON", tableName ,  "(", object@ID ,")")) )
			message(paste("Table", object@tableName, "saved as a ", object@BIO, "table") )
			}
		}
	)

setMethod("bioSave",
	signature  = "bioSaveDataFrame",
		definition = function(object) {

		tableName = paste(object@BIO, object@tableName, sep = "")

		d = object@loc

		nam = d[, object@ID]

		ranges.nam = dbGetQuery(object@CON, "select distinct bioid from ranges")$bioid

		d$has_range = is.element(nam, ranges.nam)

		res = dbWriteTable(object@CON ,tableName , d, row.names = FALSE)

		if(res) {
			dbGetQuery(object@CON,(paste("CREATE  INDEX", paste(tableName, object@ID, sep = "_") , "ON", tableName ,  "(", object@ID ,")")) )
			message(paste("Table", object@tableName, "saved as a ", object@BIO, "table") )
			} else
				message( paste("Error in saving", object@tableName) )
		}
	)

#' Import \sQuote{BIO} tables to a \code{rangeMapper} project.
#'
#' Import tables (e.g. life history data) to an active \code{rangeMapper}
#' project.
#'
#' @param con       an sqlite connection pointing to a valid \code{rangeMapper}
#'                  project.
#' @param loc       file location or \code{data.frame} name.
#' @param tableName if missing, the name of the file or \code{data.frame} is used
#' @param \dots     arguments to pass to the corresponding methods: e.g. the ID,
#'                  the column corresponding to the names of the range files
#' @return          a \sQuote{BIO} table is created in the corresponding
#'                  \code{rangeMapper} project.
#' @export          bio.save bio.merge metadata2bio
#' @examples
#'
#' require(rangeMapper)
#' require(rgdal)
#' wd = setwd(tempdir())
#' r = readOGR(system.file(package = "rangeMapper",
#' 	"extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)
#' dbcon = rangeMap.start(file = "wrens.sqlite", overwrite = TRUE,
#' 	dir = tempdir() )
#' global.bbox.save(con = dbcon, bbox = r)
#' gridSize.save(dbcon, gridSize = 2)
#' canvas.save(dbcon)
#' processRanges(spdf = r, con =  dbcon, ID = "sci_name" )
#'
#' # Upload BIO tables
#' data(wrens)
#' Troglodytes  = wrens[grep("Troglodytes", wrens$sci_name), c(2, 5)]
#' bio.save(con = dbcon, loc = Troglodytes,  ID = "sci_name")
#'
#' setwd(wd)
#'
#'
#' \dontrun{
#' require(rangeMapper)
#' require(rgdal)
#' wd = setwd(tempdir())
#' r = readOGR(system.file(package = "rangeMapper",
#'   "extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)
#' dbcon = rangeMap.start(file = "wrens.sqlite", overwrite = TRUE,
#' 	dir = tempdir() )
#' global.bbox.save(con = dbcon, bbox = r)
#' gridSize.save(dbcon, gridSize = 2)
#' canvas.save(dbcon)
#' processRanges(spdf = r, con =  dbcon, ID = "sci_name", metadata = rangeTraits() )
#'
#' wrensPath = system.file(package = "rangeMapper", "data", "wrens.csv")
#' bio.save(con = dbcon, loc = wrensPath,  ID = "sci_name")
#' bio.merge(dbcon, "wrensNew")
#' metadata2bio(dbcon)
#' setwd(wd)
#'
#' }
#'
#'
bio.save <- function(con, loc, tableName, ...) {
	if(is.character(loc)) {
		if(missing(tableName)) tableName = gsub("\\.", "_", basename(loc))
		dat = new("bioSaveFile", CON = con, loc = loc, tableName = tableName, ...) }

	if(is.data.frame(loc)) {
		if(missing(tableName)) tableName = deparse(substitute(loc))
		dat = new("bioSaveDataFrame", CON = con, loc = loc, tableName = tableName, ...) }

	bioSave(dat)

	}

#' @rdname bio.save
bio.merge <-  function(con, tableName, ...) {
	# merge 2 or more BIO tables, default is merge all

	r = new("rangeMap", CON = con)
	tableName = paste(r@BIO, tableName, sep = "")
	dots = list(...)

	if(length(dots) > 0)
	 btabs = paste(r@BIO, dots, sep = "") else
	 btabs = dbGetQuery(con, paste("select name from sqlite_master where type = 'table' and tbl_name like '", r@BIO,"%'", sep = ""))$name

	ok = sapply(btabs, function(x) dbtable.exists(con, x) )

	if(!all(ok))
		stop(paste( dQuote(names(ok[!ok])), "is not a table of this rangeMapper project"))

	ids = sapply(btabs, function(x) extract.indexed(con, x) )

	head = paste("(", paste(paste("SELECT DISTINCT",ids,"as", r@BIOID , "FROM",  names(ids)), collapse = " UNION "), ") as x")

	colnames = unlist(lapply(btabs,
		function(x) {
			a = dbGetQuery(con, paste("pragma table_info(" , shQuote(x),")" ))$name
			alias = paste(a, gsub(r@BIO, "", x) , sep = "_" )
			nm = paste(x, a,sep = "." )
			paste(nm, alias, sep = " as ")
			} ))

	colnames = paste(setdiff(colnames, ids), collapse = ",")
	colnames = paste("SELECT DISTINCT", paste(paste("x", r@BIOID, sep = "."), "as", r@BIOID ),  ",", colnames)


	joins = character()
		for(i in 1:length(ids)) joins[i] = paste("LEFT JOIN", names(ids[i]), "ON", paste("x", r@BIOID, sep ="."), "=", paste(names(ids[i]), ids[i], sep =".") )

	sqls = paste("create table", tableName,  "as", colnames, "FROM", head,  paste(joins, collapse = " " ))
	# make table
	res = dbGetQuery(con, sqls)

	# add index
	dbGetQuery(r@CON,( paste("CREATE  INDEX", paste(tableName, r@BIOID, sep = "_") , "ON", tableName ,  "(", r@BIOID ,")") ) )
    }

#' @rdname bio.save
metadata2bio <-function(con, ...) {

	r = new("rangeMap", CON = con)

	dat = dbGetQuery(r@CON, paste("select * from",  r@METADATA_RANGES) )

	if(nrow(dat) == 0) stop( paste("Empty", r@METADATA_RANGES, "table"))


	b = new("bioSaveDataFrame", CON = con, loc = dat, tableName = r@METADATA_RANGES, ID = r@BIOID, ...)

	bioSave(b)
	}

