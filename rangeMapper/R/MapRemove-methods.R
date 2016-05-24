setGeneric("rangeMapRemove", function(object, ...)   		     	standardGeneric("rangeMapRemove") )


setMethod("rangeMapRemove",
	signature = "rangeMapRemove",
	definition = function(object){

	if(length(object@tableName) == 0 )
		object@tableName = dbGetQuery(object@CON,
			'select name from sqlite_master where type = "table" and
			(tbl_name like "MAP_%" OR tbl_name like "BIO_%")')$name

	if( length(object@tablePrefix) > 0 )
	object@tableName = object@tableName[grep(object@tablePrefix, object@tableName)]

	if(length(object@tableName) > 0) {

	sql = paste("DROP TABLE IF EXISTS", object@tableName )

		for (i in 1:length(sql)) {
	  message("SQLITE:", sql[i])
	  dbGetQuery(object@CON , sql[i])
		  }
	   }
	})

#' Remove tables from a give project
#'
#' Remove tables given prefix attribute or by name
#'
#' @param  con   A valid sqlite connection.
#' @param  \dots Arguments passed to the corresponding methods specifically
#'               \sQuote{tablePrefix} or \sQuote{tableName}
#' @note         The default \sQuote{rm.rangeMapper(con)} will remove all \sQuote{MAP}
#'               and \sQuote{BIO} tables.
#' @export
#'
rm.rangeMapper <- function(con, ...) {
	 x =  new("rangeMapRemove", CON = con, ...)
	 rangeMapRemove(x)
 	}































