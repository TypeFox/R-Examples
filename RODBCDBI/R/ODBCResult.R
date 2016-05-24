#' Class ODBCResult.
#'
#' ODBC's query results class. This classes encapsulates the result of an SQL statement (either select or not). The main generator is dbSendQuery.
#' @keywords internal
#' @export
setClass(
  "ODBCResult",
  contains = "DBIResult",
  slots= list(
    connection="ODBCConnection",
    sql="character",
    state="environment"
  )
)

is_done <- function(x) {
  x@state$is_done
}
`is_done<-` <- function(x, value) {
  x@state$is_done <- value
  x
}

#' @inheritParams DBI::rownamesToColumn
#' @export
#' @rdname odbc-query
setMethod("dbFetch", "ODBCResult", function(res, n = -1, ...) {
  result <- sqlQuery(res@connection@odbc, res@sql, max=ifelse(n==-1, 0, n))
  is_done(res) <- TRUE
  result
})

#' @rdname odbc-query
#' @export
setMethod("dbHasCompleted", "ODBCResult", function(res, ...) {
  is_done(res)
})

#' @rdname odbc-query
#' @export
setMethod("dbClearResult", "ODBCResult", function(res, ...) {
  name <- deparse(substitute(res))
  is_done(res) <- FALSE
  TRUE
})


#' Database interface meta-data.
#' 
#' See documentation of generics for more details.
#' 
#' @param dbObj An object inheriting from \code{\linkS4class{ODBCConnection}}, \code{\linkS4class{ODBCDriver}}, or a \code{\linkS4class{ODBCResult}}
#' @param res An object of class \code{\linkS4class{ODBCResult}}
#' @param ... Ignored. Needed for compatibility with generic
#' @examples
#' \dontrun{
#' library(DBI)
#' data(USArrests)
#' con <- dbConnect(RODBCDBI::ODBC(), dsn="test", user="sa", password="Password12!")
#' dbWriteTable(con, "t1", USArrests, overwrite=TRUE)
#' dbWriteTable(con, "t2", USArrests, overwrite=TRUE)
#' 
#' dbListTables(con)
#' 
#' rs <- dbSendQuery(con, "select * from t1 where UrbanPop >= 80")
#' dbGetStatement(rs)
#' dbHasCompleted(rs)
#' 
#' info <- dbGetInfo(rs)
#' names(info)
#' info$fields
#' 
#' dbFetch(rs, n=2)
#' dbHasCompleted(rs)
#' info <- dbGetInfo(rs)
#' info$fields
#' dbClearResult(rs)
#' 
#' # DBIConnection info
#' names(dbGetInfo(con))
#' 
#' dbDisconnect(con)
#' }
#' @name odbc-meta
NULL

#' @rdname odbc-meta
#' @export
setMethod("dbGetRowCount", "ODBCResult", function(res, ...) {
  df <- sqlQuery(res@connection@odbc, res@sql)
  nrow(df)
})

#' @rdname odbc-meta
#' @export
setMethod("dbGetStatement", "ODBCResult", function(res, ...) {
  res@sql
})

#' @rdname odbc-meta
#' @export
setMethod("dbGetInfo", "ODBCResult", function(dbObj, ...) {
  dbGetInfo(dbObj@connection)
})


#' @rdname odbc-meta
#' @export
setMethod("dbColumnInfo", "ODBCResult", function(res, ...) {
  df <- sqlQuery(res@connection@odbc, res@sql, max=1)
  data_type <- sapply(df, class)
  data.frame(
    name=colnames(df),
    data.type=data_type,
    field.type=-1, #Can implement it(Data type in DBMS) through RODBC
    len=-1,
    precision=-1, 
    scale=-1,
    nullOK=sapply(df, function(x){any(is.null(x))}) #adhoc...
  )
})
