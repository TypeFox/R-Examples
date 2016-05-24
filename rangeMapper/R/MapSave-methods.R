
setGeneric("rangeMapSave", function(object,FUN,formula, ...)  standardGeneric("rangeMapSave") )

# method for species richness,
setMethod("rangeMapSave",
	signature = c(object = "rangeMapSave", FUN = "missing", formula = "missing"),
	definition = function(object, FUN, formula, ...){

		if(length(object@tableName) == 0) object@tableName = "species_richness"

		#build tableName
		tableName = paste(object@MAP, object@tableName, sep = "")

		# build sql subset
		sset = subsetSQLstring(object@CON, object@subset)

		# build sql string
		richnessSQL = paste("SELECT id, count(r.id) as", object@tableName ,"from ranges as r",
						if(!is.null(sset)) paste("WHERE", sset), "group by r.id")

		# build table and index
		dbGetQuery(object@CON, paste("CREATE TABLE" ,tableName, "(", object@ID, "INTEGER,",object@tableName, "NUMERIC) -- DDL:", richnessSQL))
		dbGetQuery(object@CON,paste("CREATE  INDEX", paste(object@tableName, object@ID, sep = "_") ,"ON", tableName, "(id)") )
		dbGetQuery(object@CON, paste("INSERT INTO" ,tableName, richnessSQL) )

	 	return(dbtable.exists(object@CON, tableName))


		})

# aggregate method using sqlite
setMethod("rangeMapSave",
	signature  = c(object = "rangeMapSave", FUN = "character", formula = "missing"),
		definition = function(object, FUN, formula, ...) {

		# CHECKS
		biotab = paste(object@BIO, object@biotab, sep = "")
			if(!dbtable.exists(object@CON,biotab) )
			stop( paste(sQuote(object@biotab), "is not a table of", sQuote(dbGetInfo(object@CON)$dbname)))
		# object@biotrait should exist as a field in biotab
		if(!dbfield.exists(object@CON,biotab, object@biotrait) )
			stop(paste(sQuote(object@biotrait), "is not a field of", sQuote(object@biotab)))

		initExtension(object@CON)

		# fun should  be known by sqlite
		sqlAggregate(FUN)

		# BIO_tab name
		biotab = paste(object@BIO, object@biotab, sep = "")

		#build MAP_ tableName
		tableName = paste(object@MAP, object@tableName, sep = "")

		# build sql subset
		sset = subsetSQLstring(object@CON, object@subset)
		# build sql string
		sql = paste("SELECT r.id, b.",object@biotrait,"FROM ranges r left join ",
				biotab, " b WHERE r.bioid = b.", extract.indexed(object@CON, biotab),
				  if(!is.null(sset)) paste("AND", sset) )
		sql = paste("SELECT id,", FUN ,"(", object@biotrait, ") as", object@biotrait, "from (",sql,") group by id")

		# build table and index
		dbGetQuery(object@CON, paste("CREATE TABLE",tableName, "(", object@ID, "INTEGER,",object@biotrait, "NUMERIC)"))
		dbGetQuery(object@CON, paste("CREATE INDEX", paste(tableName, "id", sep = "_") , "ON", tableName, "(id)") )
		dbGetQuery(object@CON, paste("INSERT INTO" ,tableName, sql) )


		return(dbtable.exists(object@CON, tableName))


		cat(strwrap(sql, width = 100))
		}
	)

# aggregate method using R functions called directly on the data
setMethod("rangeMapSave",
	signature  = c(object = "rangeMapSave", FUN = "function", formula = "missing"),
	definition = function(object, FUN, formula, cl, ...) {

		# tableName
		tableName = paste(object@MAP, object@tableName, sep = "")

		# get data afer checking
		dl = .rangeMapSaveData (object)
		if (missing(cl)){
			# apply R function
			X = sapply(dl, FUN = function(x) FUN(x[, object@biotrait], ...) )
		}else{
			if (clNumeric<- is.numeric(cl)) cl<- makeCluster(cl)
			# apply R function
			X = parSapply(cl=cl, dl, FUN = function(x) FUN(x[, object@biotrait], ...) )
			if (clNumeric) stopCluster(cl)
		}

		X = data.frame(id = names(X), X)
		names(X) = c(object@ID, object@biotrait)
		row.names(X) = NULL

		# build table and index
		dbGetQuery(object@CON, paste("CREATE TABLE" ,tableName, "(", object@ID, "INTEGER,",object@biotrait, "NUMERIC)"))
		dbGetQuery(object@CON, paste("CREATE INDEX", paste(tableName, "id", sep = "_") , "ON", tableName, "(id)") )
		dbWriteTable(object@CON, tableName, X, row.names = FALSE, append = TRUE)

		return(dbtable.exists(object@CON, tableName))

		})

# aggregate method using R functions called directly using formula, data interface
setMethod("rangeMapSave",
	signature  = c(object = "rangeMapSave", FUN = "function", formula = "formula"),
	definition = function(object, FUN, formula, cl, ...) {

		# tableName
		tableName = paste(object@MAP, object@tableName, sep = "")

		# get data afer checking
		dl = .rangeMapSaveData (object)

		if (missing(cl)){
		  # apply R function
		  X = sapply(dl, FUN = function(x) FUN(formula = formula, data = x, ...) )
		}else{
		  if (clNumeric<- is.numeric(cl)) cl<- makeCluster(cl)
		  # apply R function
		  X = parSapply(cl=cl, dl, FUN = function(x) FUN(formula = formula, data = x, ...) )
		  if (clNumeric) stopCluster(cl)
		}

		X = data.frame(id = names(X), X)
		names(X) = c(object@ID, object@biotrait)
		row.names(X) = NULL

		# build table and index
		dbGetQuery(object@CON, paste("CREATE TABLE" ,tableName, "(", object@ID, "INTEGER,",object@biotrait, "NUMERIC)"))
		dbGetQuery(object@CON, paste("CREATE INDEX", paste(tableName, "id", sep = "_") , "ON", tableName, "(id)") )
		dbWriteTable(object@CON, tableName, X, row.names = FALSE, append = TRUE)

		dbtable.exists(object@CON, tableName)

		})

# IMPORT RASTER
setGeneric("rangeMapImport", function(object,FUN, ...)  	  standardGeneric("rangeMapImport") )

# method for  importing external files
setMethod("rangeMapImport",
	signature  = c(object = "MapImport", FUN = "function"),
		definition = function(object,FUN, ...) {

	filenam = basename(object@path)

	if(length(object@tableName)== 0) tableName = RSQLite::make.db.names(filenam)

	tableName = paste(object@MAP, object@tableName, sep = "")

	cnv = canvas.fetch(object@CON)
	message("Converting canvas to polygons...")
	cnv = raster::rasterToPolygons(raster(cnv))

	message("Loading external MAP data...")
	rst = raster::raster(object@path)

	p4s_is_ok = proj4string_is_identical(CRSargs(CRS(proj4string(cnv))), CRSargs(raster::projection(rst, FALSE)))
	if(! p4s_is_ok )
		warning(sQuote(filenam), " may have a different PROJ4 string;\n", "canvas:", CRSargs(CRS(proj4string(cnv))), "\n", filenam, ":", CRSargs(projection(rst, FALSE)) )

	rstp = as(as(rst, "SpatialGridDataFrame"), "SpatialPointsDataFrame")
	message("Extracting Layer 1...")
	rstp = rstp[which(!is.na(rstp@data[,1])), ]

	rstp@data$ptid = as.numeric(rownames(rstp@data)) # add point id

	message(paste("Performing overlay: canvas polygons over", filenam, "...") )
	o = over(rstp, cnv)
	o$ptid = as.numeric(rownames(o))

	o = merge(o, rstp@data, by = "ptid")
	o$ptid = NULL

	message("Agregating data...")
	o = aggregate(o[, 2], list(o[,1]), FUN = FUN, na.rm = TRUE, ...)

	names(o) = c(object@ID, object@tableName)

	# build table and index
	message("Creating table and indexes...")
	dbGetQuery(object@CON, paste("CREATE TABLE" ,tableName, "(", object@ID, "INTEGER,",object@tableName, "FLOAT)"))
	dbGetQuery(object@CON, paste("CREATE INDEX", paste(tableName, "id", sep = "_") , "ON", tableName, "(id)") )
	dbWriteTable(object@CON, tableName, o, row.names = FALSE, append = TRUE)


	res = dbtable.exists(object@CON, tableName)

	if(res) message(paste(sQuote(basename(object@path)), "imported"))

	return(res)


		}
	)

#' Save, retrieve and export maps.
#'
#' Apply a chosen \code{SQL} or function at each grid cell, allowing for
#' complex subsetting at both ID (e.g. species) and pixel (e.g. assemblage)
#' levels.
#'
#' The subset argument accepts a named list. Names refers to \sQuote{BIO},
#' \sQuote{MAP} and \sQuote{metadata_rages} table names while the strings in
#' the list are character strings containing the SQL \code{WHERE} clause. The
#' subset can point to either one table type (e.g.
#' \code{list(MAP_species_richness = "species_richness > 500")} ) or can point
#' to several table types (e.g. \code{list(BIO_lifeHistory = "clutch_size > 4",
#' MAP_meanAltitude = "meanAltitude < 1000", metadata_ranges = "Area < 1000")})
#'
#' Any valid SQL expression can be used to build up a subset. See
#' \url{http://www.sqlite.org/lang_expr.html}
#'
#' When using \code{cl} parameter you must load the apropiated packages used in
#' \code{FUN} by loading the packages inside the function or initializing the
#' cluster before calling rangeMap.save (e.g. \code{clusterEvalQ(cl=cl, library(caper))})).
#'
#' @param CON       an sqlite connection pointing to a valid \code{rangeMapper}
#'                  project.
#' @param tableName name of the table (quoted) to be added to the sqlite database.
#'                  the prefix \sQuote{MAP} will be appended to \code{tableName}
#'                  prior to saving.
#' @param FUN       the function to be applied to each pixel. If \code{FUN} is
#'                  missing then species richness (species count) is computed.
#' @param biotab    character string identifying the \sQuote{BIO} table to use.
#' @param biotrait  character string identifying the ID of the \sQuote{BIO}
#'                  table. see \code{\link{bio.save}}
#' @param subset    a named \code{\link{list}}. See details
#' @param path      path to the raster file(quoted) to be imported to the existing
#'                  project. \code{raster package} is required at this step.
#' @param overwrite if \code{TRUE} then the table is removed
#' @param cl        the number of cores to use or a cluster object defined with
#'                  \code{\link[parallel]{makeCluster}} in package \code{parallel}
#'                  or \code{\link[snow]{makeCluster}} from \code{snow} package.
#' @param \dots     when \code{FUN} is a function, \dots{} denotes any extra
#'                  arguments to be passed to it.
#' @return          \code{TRUE} when the MAP was created successfully.
#' @note            \code{SQL} aggregate functions are more efficient then their R
#'                  counterparts. For simple aggregate functions like mean, median, sd, count
#'                  it is advisable to use \code{SQL} functions rather then R functions.
#' @seealso         \code{\link{metadata.update}}.
#' @export
#' @examples
#' require(rangeMapper)
#' require(rgdal)
#' breding_ranges = readOGR(system.file(package = "rangeMapper",
#'      "extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)
#' breding_ranges = spTransform(breding_ranges,
#'     CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0
#'         +ellps=WGS84 +units=m +no_defs") )
#' data(wrens)
#' d = subset(wrens, select = c('sci_name', 'body_size', 'body_mass', 'clutch_size') )
#'
#' con = ramp("wrens.sqlite", gridSize = 200000, spdf = breding_ranges, biotab = d,
#'             ID = "sci_name", metadata = rangeTraits(),
#'             FUN = "median", overwrite = TRUE)
#'
#'
#' lmSlope = function(formula, data) {
#'     fm = try(lm(formula, data = data), silent = TRUE)
#'     if (inherits(fm, "try-error"))
#'         res = NA else res = coef(fm)[2]
#'     as.numeric(res)
#' }
#'
#' # Subsetting by Species and Assembladge
#' rangeMap.save(con, FUN = lmSlope, biotab = "biotab", biotrait = "body_mass",
#'     tableName = "slope_bodyMass_clutchSize", formula = log(body_mass) ~ clutch_size,
#'     list(MAP_species_richness = "species_richness >= 5",
#'         BIO_biotab = "body_size > 15"
#'         ), overwrite = TRUE)
#'
#' \dontrun{
#' # Import raster maps to the current project
#' r = system.file(package = "rangeMapper", "extdata", "etopo1", "etopo1_Americas.tif")
#' rangeMap.save(con, path = r, tableName = "meanAltitude", FUN = mean, overwrite = TRUE)
#' }
#'
#'
#' m = rangeMap.fetch(con, spatial = FALSE)
#' plot(m)
#'
rangeMap.save  <- function(CON, tableName, FUN, biotab, biotrait, subset = list(), path , overwrite = FALSE, cl, ...) {

	if(overwrite & !missing(tableName))
	try(dbGetQuery(CON, paste("DROP TABLE", paste("MAP", tableName, sep = "_"))), silent = TRUE)

	if(overwrite & missing(tableName))
	try(dbGetQuery(CON, "DROP TABLE MAP_species_richness"), silent = TRUE)



	o = FALSE

	if(!missing(path)) { #  external map
			if(missing(tableName))
				rmap = new("MapImport", CON = CON, path = path) else
				rmap = new("MapImport", CON = CON, path = path, tableName = tableName)

			 o = rangeMapImport(rmap, FUN = FUN)
			}

	if(missing(FUN) ) { #species richness
			if(missing(tableName))
				rmap = new("rangeMapSave", CON = CON, subset = subset) else
				rmap = new("rangeMapSave", CON = CON, tableName = tableName, subset = subset)
			 o = rangeMapSave(rmap)

			}

	if(!missing(FUN) & missing(path)) { # SQL or R function
			rmap = new("rangeMapSave", CON = CON,  biotab = biotab, biotrait = biotrait, tableName = tableName, subset = subset)
			 o = rangeMapSave(rmap, FUN = FUN, ...)
			}
		o
		}


