#' ODBCDriver and methods.
#' 
#' An ODBC driver implementing the R database (DBI) API.
#' This class should always be initialized with the \code{ODBC()} function.
#' It returns an object that allows you to connect to ODBC.
#'
#' @export
#' @keywords internal
setClass("ODBCDriver", contains = "DBIDriver")

#' Generate an object of ODBCDriver class
#'
#' This driver is for implementing the R database (DBI) API.
#' This class should always be initialized with the \code{ODBC()} function.
#' ODBC driver does nothing for ODBC connection. It just exists for S4 class compatibility with DBI package. 
#'
#' @export
#' @examples
#' \dontrun{
#' driver <- RODBCDBI::ODBC()
#' # Connect to a ODBC data source
#' con <- dbConnect(driver, dsn="test")
#' # Always cleanup by disconnecting the database
#' #' dbDisconnect(con)
#' }
ODBC <- function() {new("ODBCDriver")}


#' @rdname ODBCDriver-class
#' @export
setMethod("dbUnloadDriver", "ODBCDriver", function(drv, ...) {TRUE})

#' Connect/disconnect to a ODBC data source
#'
#' These methods are straight-forward implementations of the corresponding generic functions.
#' 
#' @param drv an object of class ODBCDriver
#' @param dsn Data source name you defined by ODBC data source administrator tool.
#' @param user User name to connect as.
#' @param password Password to be used if the DSN demands password authentication.
#' @param ... Other parameters passed on to methods
#' @export
#' @examples
#' \dontrun{
#' # Connect to a ODBC data source
#' con <- dbConnect(RODBCDBI::ODBC(), dsn="test")
#' # Always cleanup by disconnecting the database
#' #' dbDisconnect(con)
#' }
setMethod(
  "dbConnect", 
  "ODBCDriver", 
  function(drv, dsn, user = NULL, password = NULL, ...){
    uid <- if(is.null(user)) "" else user
    pwd <- if(is.null(password)) "" else password
    connection <- odbcConnect(dsn, uid, pwd, ...)
    new("ODBCConnection", odbc=connection)
  }
)

#' @rdname ODBCDriver-class
#' @export
setMethod("dbIsValid", "ODBCDriver", function(dbObj) {TRUE})

#' Get ODBCDriver metadata.
#' 
#' Nothing to do for ODBCDriver case
#' 
#' @rdname ODBCDriver-class
#' @export
setMethod("dbGetInfo", "ODBCDriver", function(dbObj, ...) {NULL})
