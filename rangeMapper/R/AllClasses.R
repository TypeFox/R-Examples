
setClass("rangeMapStart",
	representation(
		dir = "character",
		file = "character",
		skeleton = "list",
		overwrite = "logical"
		),
	prototype(
		file = paste("rangeMapperProj",format(Sys.time(), "%Y-%m-%d_%H-%M-%S.sqlite"), sep = "_"),
		skeleton = 	list(
	create =
	c(	"CREATE TABLE version (ver CHAR)",
		"CREATE TABLE proj4string (p4s CHAR)",
		"CREATE TABLE gridSize(gridSize FLOAT)",
		"CREATE TABLE bbox(xmin FLOAT,xmax FLOAT,ymin FLOAT,ymax FLOAT)",
		"CREATE TABLE canvas (x FLOAT,y FLOAT,id INT)",
		"CREATE TABLE ranges (id INT,bioid CHAR)",
		"CREATE TABLE metadata_ranges (bioid CHAR)"

	),
	index =
	c(  "CREATE INDEX IF NOT EXISTS   id_canvas ON canvas (id)",
		"CREATE INDEX IF NOT EXISTS   id_ranges ON ranges (id)",
		"CREATE INDEX IF NOT EXISTS   bioid_ranges ON ranges (bioid)",
		"CREATE INDEX IF NOT EXISTS   bioid_metadata_ranges ON metadata_ranges (bioid)")
	),
	overwrite = FALSE
	),

	validity = function(object) {
	if(!file.exists(object@dir)) stop("'dir' should be set and point to a valid location.")
	} )

setClass("rangeMap",
	representation(
		CON = "SQLiteConnection",
		VERSION = "character",
		ID = "character",
		BIOID = "character",
		PROJ4STRING = "character",
		GRIDSIZE = "character",
		BBOX = "character",
		METADATA_RANGES = "character",
		CANVAS = "character",
		RANGES = "character",
		BIO = "character",
		MAP = "character"
		),
	prototype(
		VERSION = "version",
		ID = "id",
		BIOID = "bioid",
		PROJ4STRING = "proj4string",
		GRIDSIZE = "gridsize",
		BBOX = "bbox",
		METADATA_RANGES = "metadata_ranges",
		CANVAS = "canvas",
		RANGES = "ranges",
		BIO = "BIO_",
		MAP = "MAP_"
		),

	validity = function(object) {

	if ( ! all( dbtable.exists(object@CON, object@PROJ4STRING),
			  dbtable.exists(object@CON, object@METADATA_RANGES),
			  dbtable.exists(object@CON, object@GRIDSIZE),
			  dbtable.exists(object@CON, object@BBOX),
			  dbtable.exists(object@CON, object@CANVAS),
			  dbtable.exists(object@CON, object@RANGES) ) ) stop ("Corrupt rangeMapper project!")



	})

setClass("rangeFiles",
	representation(
	dir = "character",
	ogr = "logical",
	polygons.only = "logical") ,
	prototype(
	dir = "",
	ogr = TRUE,
	polygons.only = TRUE),
	validity = function(object) {
	if(!file.exists(object@dir)) stop("'dir' should be set and point to a valid location.")

	})

setClass("gridSize",
	representation(
	gridSize = "numeric"
		),

	contains = "rangeMap",

	validity = function(object)	{
		if(!is.empty(object@CON, object@GRIDSIZE)) stop("The grid size was already set!")
		if(is.empty(object@CON, object@BBOX)) stop("There is no bounding box!")

	})

setClass("rangeMapSave",
	representation(
		biotab    = "character",
		biotrait  = "character",
		tableName = "character",
		subset    = "list"),

	contains = "rangeMap",

	validity = function(object)	{
	# the new table should not exist
		if(dbtable.exists(object@CON, paste(object@MAP, object@tableName, sep = "") ) )
			stop(paste(sQuote(object@tableName), " already exists."))

		if(dbtable.exists(object@CON, paste(object@BIO, object@tableName, sep = "") ) )
			stop( paste(sQuote(object@tableName), " already exists."))

	})

setClass("bioSaveFile",
	representation(loc = "character", sep = "character"),
							contains = "rangeMapSave",
	prototype( sep = ";", tableName = "unknown"),
	validity = function(object) {
					if(!file.exists(object@loc))
					stop(paste(sQuote(object@loc), "is not a valid file"))
							})

setClass("bioSaveDataFrame", representation(loc = "data.frame"),
	contains = "rangeMapSave",
	validity = function(object) {
	})

setClass("MapImport", representation(path = "character"),
	contains = "rangeMapSave",

	validity = function(object) {

	if(!file.exists(object@path)) stop(paste(sQuote(object@path), "is not a valid path.") )
	if(!require("raster")) stop("package raster is not available")

	})

setClass("rangeMapFetch", representation(
	tableName    = "character"),
	contains     = "rangeMap",

	validity = function(object)	{
	mapNam = paste(object@MAP, object@tableName, sep = "")

	invalidNam = sapply(mapNam, FUN = function(x) !dbtable.exists(object@CON, x) )

	if( any(invalidNam) )
	  stop(paste(sQuote(names(invalidNam[invalidNam] )), "is not a valid MAP!\n"))

	# check if empty map
	mapvar = sapply(mapNam, function(x)
				setdiff(dbGetQuery(object@CON, paste("pragma table_info(", x, ")"))$name, object@ID ) )

	sql = paste("select count (", mapvar, ") FROM", mapNam)
	isempty = sapply(sql, function(x) dbGetQuery(object@CON,  x)[, 1] ) < 1

	if(any(isempty))
	 stop(paste(sQuote(mapNam[isempty]), "is an empty MAP!\n"))
	})

setClass("rangeMapRemove", representation(
	tableName = "character",
	tablePrefix = "character"
	),
	 contains = "rangeMap",

	validity = function(object)	{
				if(object@PROJ4STRING%in%object@tableName |
				   object@BBOX%in%object@tableName|
				   object@GRIDSIZE%in%object@tableName|
				   object@CANVAS%in%object@tableName|
				   object@CANVAS%in%object@tableName)
	stop( paste(object@PROJ4STRING,",", object@BBOX, ",",object@GRIDSIZE,",",
					object@CANVAS, "or",
					object@RANGES, "table(s) cannot be removed!"))
	})

setClass("SpatialPixelsRangeMap", representation(
	mapvar    = "character"),
    contains  = "SpatialPixelsDataFrame",

	validity = function(object)	{
	return(TRUE)
	})

setOldClass(c('data.table'))
setOldClass(c('data.table', 'rmap.frame'))
setAs("data.table", "rmap.frame", function(from) {
    as.rmap.frame(from)
	})