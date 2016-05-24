#' Send query to SQLite database
#' 
#' @param db_name character. Path to database.
#' @param query_string character. Query string.
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbSendQuery dbGetQuery dbWriteTable
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
sendQuery = function(db_name, query_string) {
	db = RSQLite::dbConnect(RSQLite::SQLite(), db_name)
	tryCatch({
				RSQLite::dbSendQuery(db, query_string)
			}, finally = {
				RSQLite::dbDisconnect(db)
			})
}

#' Get query results from a SQLite database
#' 
#' @param db_name character. Path to database.
#' @param query_string character. Query string.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
getQuery = function(db_name, query_string) {
	db = RSQLite::dbConnect(RSQLite::SQLite(), db_name)
	tryCatch({
				return(RSQLite::dbGetQuery(db, query_string))
			}, finally = {
				RSQLite::dbDisconnect(db)
			})
}