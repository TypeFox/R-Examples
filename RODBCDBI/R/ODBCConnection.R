#' @include RODBCDBI.R
NULL

#' Class ODBCConnection.
#'
#' \code{ODBCConnection} objects are usually created by \code{\link[DBI]{dbConnect}}
#' @keywords internal
#' @export
setClass(
  "ODBCConnection",
  contains="DBIConnection",
  slots=list(odbc="RODBC")
)

#' Execute a statement on a given database connection.
#' 
#' To retrieve results a chunk at a time, use \code{dbSendQuery},
#' \code{dbFetch}, then \code{ClearResult}. Alternatively, if you want all the
#' results (and they'll fit in memory) use \code{dbGetQuery} which sends,
#' fetches and clears for you.
#' 
#' @param conn An existing \code{\linkS4class{ODBCConnection}}
#' @param statement  The SQL which you want to run
#' @param res An object of class \code{\linkS4class{ODBCResult}}
#' @param n Number of rows to return. If less than zero returns all rows.
#' @param ... Other parameters passed on to methods
#' @export
#' @rdname odbc-query
setMethod("dbSendQuery", "ODBCConnection", function(conn, statement, ...) {
  statement <- enc2utf8(statement)  
  env <- new.env()
  assign("is_done", FALSE, envir=env)
  new("ODBCResult", connection=conn, sql=statement, state=env)
})

#' Get DBMS metadata.
#' 
#' @param dbObj An object inheriting from \code{\linkS4class{ODBCConnection}}, \code{\linkS4class{ODBCDriver}}, or a \code{\linkS4class{ODBCResult}}
#' @param ... Other parameters passed on to methods
#' @export
setMethod("dbGetInfo", "ODBCConnection", function(dbObj, ...) {
  info <- RODBC::odbcGetInfo(dbObj@odbc)
  list(dbname = unname(info["DBMS_Name"]), 
       db.version = unname(info["DBMS_Ver"]), 
       username = "", 
       host = "", 
       port = "",
       sourcename = unname(info["Data_Source_Name"]),
       servername = unname(info["Server_Name"]),
       drivername = unname(info["Driver_Name"]),
       odbc.version = unname(info["ODBC_Ver"]),
       driver.version = unname(info["Driver_Ver"]),
       odbcderiver.version = unname(info["Driver_ODBC_Ver"]))
})


#' List fields in specified table.
#' 
#' @param conn An existing \code{\linkS4class{ODBCConnection}}
#' @param name a length 1 character vector giving the name of a table.
#' @export
#' @examples
#' \dontrun{
#' library(DBI)
#' con <- dbConnect(RODBCDBI::ODBC(), dsn="test", user="sa", password="Password12!")
#' dbWriteTable(con, "iris", iris, overwrite=TRUE)
#' dbListFields(con, "iris")
#' dbDisconnect(con)
#' }
setMethod("dbListFields", c("ODBCConnection", "character"), function(conn, name) {
  sqlColumns(conn@odbc, name)$COLUMN_NAME
})

#' List available ODBC tables.
#' 
#' @param conn An existing \code{\linkS4class{ODBCConnection}}
#' @export
setMethod("dbListTables", "ODBCConnection", function(conn){
  sqlTables(conn@odbc)$TABLE_NAME
})

#' Write a local data frame or file to the database.
#' 
#' @export
#' @rdname dbWriteTable
#' @param conn a \code{\linkS4class{ODBCConnection}} object, produced by \code{\link[DBI]{dbConnect}}
#' @param name a character string specifying a table name. ODBCConnection table names 
#'   are \emph{not} case sensitive, e.g., table names \code{ABC} and \code{abc} 
#'   are considered equal.
#' @param value a data.frame (or coercible to data.frame) object or a 
#'   file name (character).  when \code{value} is a character, it is interpreted as a file name and its contents imported to ODBC.
#' @param overwrite logical. Should data be overwritten?
#' @param append logical. Should data be appended to an existing table?
#' @param ... additional arguments passed to the generic.
#' @export
#' @examples
#' \dontrun{
#' library(DBI)
#' con <- dbConnect(RODBCDBI::ODBC(), dsn="test", user="sa", password="Password12!")
#' dbWriteTable(con, "mtcars", mtcars, overwrite=TRUE)
#' dbReadTable(con, "mtcars") 
#' dbDisconnect(con)
#' }
setMethod("dbWriteTable", c("ODBCConnection", "character", "data.frame"), function(conn, name, value, overwrite=FALSE, append=FALSE, ...) {
  sqlSave(conn@odbc, dat=value, tablename=name, safer=!overwrite, append=append, ...)
  invisible(TRUE)
})

#' Does the table exist?
#' 
#' @param conn An existing \code{\linkS4class{ODBCConnection}}
#' @param name String, name of table. Match is case insensitive.
#' @return boolean value which indicated whether the table exists or not
#' @export
setMethod("dbExistsTable", c("ODBCConnection", "character"), function(conn, name) {
  tolower(name) %in% tolower(dbListTables(conn))
})

#' Remove a table from the database.
#' 
#' Executes the SQL \code{DROP TABLE}.
#' 
#' @param conn An existing \code{\linkS4class{ODBCConnection}}
#' @param name character vector of length 1 giving name of table to remove
#' @export
setMethod("dbRemoveTable", c("ODBCConnection", "character"), function(conn, name) {
  if(dbExistsTable(conn, name)){
    sqlDrop(conn@odbc, name)
  }
  invisible(TRUE)
})

#' Convenience functions for importing/exporting DBMS tables
#' 
#' These functions mimic their R/S-Plus counterpart \code{get}, \code{assign},
#' \code{exists}, \code{remove}, and \code{objects}, except that they generate
#' code that gets remotely executed in a database engine.
#' 
#' @return A data.frame in the case of \code{dbReadTable}; otherwise a logical
#' indicating whether the operation was successful.
#' @note Note that the data.frame returned by \code{dbReadTable} only has
#' primitive data, e.g., it does not coerce character data to factors.
#' 
#' @param conn a \code{\linkS4class{ODBCConnection}} object, produced by \code{\link[DBI]{dbConnect}}
#' @param name a character string specifying a table name.
#' @param row.names a character string specifying a table name.
#' @param check.names If \code{TRUE}, the default, column names will be converted to valid R identifiers.
#' @param select.cols  A SQL statement (in the form of a character vector of 
#'    length 1) giving the columns to select. E.g. "*" selects all columns, 
#'    "x,y,z" selects three columns named as listed.
#' @inheritParams DBI::rownamesToColumn
#' @export
#' @examples
#' \dontrun{
#' library(DBI)
#' con <- dbConnect(RODBCDBI::ODBC(), dsn="test", user="sa", password="Password12!")
#' dbWriteTable(con, "mtcars", mtcars, overwrite=TRUE)
#' dbReadTable(con, "mtcars")
#' dbGetQuery(con, "SELECT * FROM mtcars WHERE cyl = 8")
#' 
#' # Supress row names
#' dbReadTable(con, "mtcars", row.names = FALSE)
#' dbGetQuery(con, "SELECT * FROM mtcars WHERE cyl = 8", row.names = FALSE)
#' 
#' dbDisconnect(con)
#' }
setMethod("dbReadTable", c("ODBCConnection", "character"), function(conn, name, row.names = NA, check.names = TRUE, select.cols = "*") {
  out <- dbGetQuery(conn, paste("SELECT", select.cols, "FROM", name), row.names = row.names)
  if (check.names) {
    names(out) <- make.names(names(out), unique = TRUE)
  }  
  out
})

#' Close a current session.
#' 
#' @rdname dbDisconnect
#' @param conn a \code{\linkS4class{ODBCConnection}} object, produced by \code{\link[DBI]{dbConnect}}
#' @examples
#' \dontrun{
#' library(DBI)
#' con <- dbConnect(RODBCDBI::ODBC(), dsn="test", user="sa", password="Password12!")
#' dbDisconnect(con)
#' }
#' @export
setMethod("dbDisconnect", "ODBCConnection", function(conn) {
  if (RODBC:::odbcValidChannel(conn@odbc)){
    odbcClose(conn@odbc)
  } else{
    TRUE
  }
})
