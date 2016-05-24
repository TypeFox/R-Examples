#' Send query to the elliptic database
#' 
#' This API is used for write operations such as CREATE and INSERT.
#'
#' @param db an object of ecdb class
#' @param statement character, the SQL statement
#' @param ... database-specific parameters may be specified here
#'
#' @return a result set object
#'
#' @keywords ecdb
#'
#' @author Stephen H-T. Lihn
#'
#' @export 
#' 
### <======================================================================>
"ecdb.dbSendQuery" <- function(db, statement, ...)
{
    rs <- RSQLite::dbSendQuery(db@conn, statement, ...)
    RSQLite::dbClearResult(rs)
    return(rs)
}
### <---------------------------------------------------------------------->
