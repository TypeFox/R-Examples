# recursive function to check depth of list. hat tip flodel 
# at stackoverflow: http://stackoverflow.com/questions/13432863/determine-level-of-nesting-in-r
#' Determine depth of a list
#' 
#' This function recursively checks the depth of a list and returns an integer value of depth
#' 
#' @param this an object of class list 
#' @details Hat tip to flodel at stackoverflow for suggesting this light weight way analyze depth of a nested list. Further complexity needs to be added to cover the fact that data.frame are lists, too. A more sophisticated recursive function can be found in the gatveys2 package.
#' @references http://stackoverflow.com/questions/13432863/determine-level-of-nesting-in-r
#' @export
getListDepth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, getListDepth)), 0L)

#' Check Validity of a PostgreSQL connection
#' 
#' Is the PostgreSQL connection expired?
#' 
#' @param dbObj PostgreSQL connection object.
#' @importFrom DBI dbIsValid
#' @importFrom DBI dbGetInfo
#' @import RPostgreSQL
#' @import methods
#' @docType methods
#' @aliases dbIsValid
#' @rdname dbIsValid-methods
#' @export
setMethod("dbIsValid", "PostgreSQLConnection", function(dbObj) {
  isValid <- tryCatch({DBI::dbGetInfo(dbObj)},
                      error = function(e) NULL)
  !is.null(isValid)  
})


