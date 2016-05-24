#' Protective commit
#' 
#' Protective commit after sending query to the elliptic database. 
#'
#' @param db an object of ecdb class
#'
#' @return The db object
#'
#' @keywords ecdb
#'
#' @author Stephen H-T. Lihn
#'
#' @export 
#' 
### <======================================================================>
ecdb.protectiveCommit <- function(db) {
    # the commit is to prevent a strange error to happen
    tryCatch({ RSQLite::dbCommit(db@conn) }, error= function(e) {"ignore"})
    invisible(db)
}
### <---------------------------------------------------------------------->
