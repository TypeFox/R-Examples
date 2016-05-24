#' Executes the specified query and returns a data frame. This function currently
#' supports RODBC, RSQLite, and RMySQL. For other databases, use getQuery() and
#' execute the SQL statement using the appropriate database connection.
#' 
#' @param query the query to execute.
#' @param connection the database connection.
#' @param maxLevels the maximum number of levels a factor can have before being
#'        converted to a character. Set to \code{NULL} to not recode.
#' @param ... other parameters passed to \code{\link{getSQL}} and \code{\link{sqlexec}}.
#' @seealso sqlexec, cacheQuery
#' @export
execQuery <- function(query=NULL, connection=NULL, maxLevels=20, ...) {
	sql = getSQL(query=query, ...)
	df <- sqlexec(connection, sql=sql, ...)
	if(!is.null(maxLevels)) {
		df <- recodeColumns(df, maxLevels)
	}
	return(df)
}
