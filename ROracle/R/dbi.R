#
# Copyright (c) 2011, 2015, Oracle and/or its affiliates. All rights reserved.
#
#    NAME
#      dbi.R - DBI implementation for Oracle RDBMS based on OCI
#
#    DESCRIPTION
#      DBI implementation for Oracle RDBMS based on OCI.
#
#    NOTES
#
#    MODIFIED   (MM/DD/YY)
#    rpingte     03/31/15 - add ora.attributes
#    rpingte     01/29/15 - add unicode_as_utf8
#    ssjaiswa    12/16/14 - [16907374] Add date argument in .oci.WriteTable
#    ssjaiswa    09/10/14 - Add bulk_write
#    rpingte     03/11/14 - add end of result
#    rkanodia    10/03/13 - Add session mode
#    rkanodia    08/30/13 - [17383542] Added schema argument to
#                           .oci.WriteTable() & .oci.RemoveTable()
#    rkanodia    12/10/12 - Changed default value of bulk_read to 1000
#    paboyoun    09/17/12 - add difftime support
#    demukhin    09/11/12 - add Extproc driver
#    jfeldhau    06/19/12 - ROracle support for TimesTen.
#    rpingte     04/25/12 - use externalptr for connection and result set
#    rpingte     04/23/12 - add interrupt enable to driver
#    demukhin    01/20/12 - cleanup
#    paboyoun    01/04/12 - minor code cleanup
#    demukhin    12/08/11 - more OraConnection and OraResult methods
#    demukhin    12/01/11 - add support for more methods
#    demukhin    10/12/11 - creation
#

##
## Class: DBIDriver
##
setClass("OraDriver",
  representation(
    handle = "externalptr"),
  contains = "DBIDriver"
)

Oracle <- function(interruptible = FALSE, unicode_as_utf8 = TRUE,
                   ora.attributes = FALSE)
{
  .oci.Driver(.oci.drv(), interruptible = interruptible,
              unicode_as_utf8 = unicode_as_utf8,
              ora.attributes = ora.attributes)
}

setMethod("dbUnloadDriver",
signature(drv = "OraDriver"),
function(drv, ...) .oci.UnloadDriver(.oci.drv())
)

setMethod("dbListConnections",
signature(drv = "OraDriver"),
function(drv, ...) .oci.DriverInfo(.oci.drv(), "connections")[[1L]]
)

setMethod("dbGetInfo",
signature(dbObj = "OraDriver"),
function(dbObj, ...) .oci.DriverInfo(.oci.drv(), ...)
)

setMethod("summary",
signature(object = "OraDriver"),
function(object, ...) .oci.DriverSummary(.oci.drv())
)

setMethod("show",
signature(object = "OraDriver"),
function (object)
{
  .oci.DriverSummary(.oci.drv())
  invisible()
}
)

##
## Class: DBIDriver
##
setClass("ExtDriver",
  representation(
    handle = "externalptr"),
  contains = "DBIDriver"
)

Extproc <- function(extproc.ctx = NULL)
{
  .oci.Driver(.ext.drv(), extproc.ctx = extproc.ctx)
}

setMethod("dbUnloadDriver",
signature(drv = "ExtDriver"),
function(drv, ...) .oci.UnloadDriver(.ext.drv())
)

setMethod("dbListConnections",
signature(drv = "ExtDriver"),
function(drv, ...) .oci.DriverInfo(.ext.drv(), "connections")[[1L]]
)

setMethod("dbGetInfo",
signature(dbObj = "ExtDriver"),
function(dbObj, ...) .oci.DriverInfo(.ext.drv(), ...)
)

setMethod("summary",
signature(object = "ExtDriver"),
function(object, ...) .oci.DriverSummary(.ext.drv())
)

setMethod("show",
signature(object = "ExtDriver"),
function (object)
{
  .oci.DriverSummary(.ext.drv())
  invisible()
}
)

##
## Class: DBIConnection
##
setClass("OraConnection",
  representation(
    handle = "externalptr",
    timesten = "logical"),
  contains = "DBIConnection"
)

setMethod("dbConnect",
signature(drv = "OraDriver"),
function(drv, username = "", password = "", dbname = "", prefetch = FALSE,
         bulk_read = 1000L, bulk_write= 1000L , stmt_cache = 0L,
         external_credentials = FALSE, sysdba = FALSE, ...)
.oci.Connect(.oci.drv(), username = username, password = password,
             dbname = dbname, prefetch = prefetch, bulk_read = bulk_read,
             bulk_write = bulk_write, stmt_cache = stmt_cache,
             external_credentials = external_credentials, sysdba = sysdba)
)

setMethod("dbConnect",
signature(drv = "ExtDriver"),
function(drv, prefetch = FALSE, bulk_read = 1000L, bulk_write= 1000L,
         stmt_cache = 0L, external_credentials = FALSE, sysdba = FALSE, ...)
.oci.Connect(.ext.drv(), prefetch = prefetch, bulk_read = bulk_read,
             bulk_write = bulk_write, stmt_cache = stmt_cache,
             external_credentials = external_credentials, sysdba = sysdba)
)

setMethod("dbDisconnect",
signature(conn = "OraConnection"),
function(conn, ...) .oci.Disconnect(conn)
)

setMethod("dbSendQuery",
signature(conn = "OraConnection", statement = "character"),
function(conn, statement, data = NULL, prefetch = FALSE, 
         bulk_read = 1000L, bulk_write = 1000L, ...)
.oci.SendQuery(conn, statement, data = data, prefetch = prefetch,
               bulk_read = bulk_read , bulk_write = bulk_write)
)

setMethod("dbGetQuery",
signature(conn = "OraConnection", statement = "character"),
function(conn, statement, data = NULL, prefetch = FALSE, 
         bulk_read = 1000L, bulk_write = 1000L, ...)
.oci.GetQuery(conn, statement, data = data, prefetch = prefetch,
              bulk_read = bulk_read, bulk_write = bulk_write)
)

setMethod("dbGetException",
signature(conn = "OraConnection"),
function(conn, ...) .oci.GetException(conn)
)

setMethod("dbListResults",
signature(conn = "OraConnection"),
function(conn, ...) .oci.ConnectionInfo(conn, "results")[[1L]]
)

setMethod("dbGetInfo",
signature(dbObj = "OraConnection"),
function(dbObj, what, ...)
{
  if (missing(what))
    .oci.ConnectionInfo(dbObj)
  else
    .oci.ConnectionInfo(dbObj, what = what)
}
)

setMethod("summary",
signature(object = "OraConnection"),
function(object, ...) .oci.ConnectionSummary(object)
)

setMethod("show",
signature(object = "OraConnection"),
function (object)
{
  .oci.ConnectionSummary(object)
  invisible()
}
)

##
## DBIConnection: Convenience methods
##
setMethod("dbListTables",
signature(conn = "OraConnection"),
function(conn, schema = NULL, all = FALSE, full = FALSE, ...)
.oci.ListTables(conn, schema = schema, all = all, full = full)
)

setMethod("dbReadTable",
signature(conn = "OraConnection", name = "character"),
function(conn, name, schema = NULL, row.names = NULL, ...)
.oci.ReadTable(conn, name, schema = schema, row.names = row.names)
)

setMethod("dbWriteTable",
signature(conn = "OraConnection", name = "character", value = "data.frame"),
function(conn, name, value, row.names = FALSE, overwrite = FALSE,
         append = FALSE, ora.number = TRUE, schema = NULL, date = FALSE, ...)
.oci.WriteTable(conn, name, value, row.names = row.names,
                overwrite = overwrite, append = append,
                ora.number = ora.number, schema = schema, date = date)
)

setMethod("dbExistsTable",
signature(conn = "OraConnection", name = "character"),
function(conn, name, schema = NULL, ...)
.oci.ExistsTable(conn, name, schema = schema)
)

setMethod("dbRemoveTable",
signature(conn = "OraConnection", name = "character"),
function(conn, name, purge = FALSE, schema = NULL, ...)
.oci.RemoveTable(conn, name, purge = purge, schema = schema)
)

setMethod("dbListFields",
signature(conn = "OraConnection", name = "character"),
function(conn, name, schema = NULL, ...)
.oci.ListFields(conn, name, schema = schema)
)

##
## DBIConnection: Transaction management
##
setMethod("dbCommit",
signature(conn = "OraConnection"),
function(conn, ...) .oci.Commit(conn, ...)
)

setMethod("dbRollback",
signature(conn = "OraConnection"),
function(conn, ...) .oci.Rollback(conn, ...)
)

##
## DBIConnection: Stored procedures
##
setMethod("dbCallProc",
signature(conn = "OraConnection"),
function(conn, ...) .NotYetImplemented()
)

##
## Class: DBIResult
##
setClass("OraResult",
  representation(
    handle = "externalptr"),
  contains = "DBIResult"
)

setMethod("fetch",
signature(res = "OraResult"),
function(res, n = -1, ...) .oci.fetch(res, as.integer(n))
)

setMethod("dbClearResult",
signature(res = "OraResult"),
function(res, ...) .oci.ClearResult(res)
)

setMethod("dbColumnInfo",
signature(res = "OraResult"),
function(res, ...) .oci.ResultInfo(res, "fields")[[1L]]
)

setMethod("dbGetStatement",
signature(res = "OraResult"),
function(res, ...) .oci.ResultInfo(res, "statement")[[1L]]
)

setMethod("dbHasCompleted",
signature(res = "OraResult"),
function(res) .oci.EOFResult(res)
)

setMethod("dbGetRowsAffected",
signature(res = "OraResult"),
function(res, ...) .oci.ResultInfo(res, "rowsAffected")[[1L]]
)

setMethod("dbGetRowCount",
signature(res = "OraResult"),
function(res, ...) .oci.ResultInfo(res, "rowCount")[[1L]]
)

setMethod("dbGetInfo",
signature(dbObj = "OraResult"),
function(dbObj, what, ...)
{
  if (missing(what))
    .oci.ResultInfo(dbObj)
  else
    .oci.ResultInfo(dbObj, what = what)
}
)

setMethod("summary",
signature(object = "OraResult"),
function(object, ...) .oci.ResultSummary(object)
)

setMethod("show",
signature(object = "OraResult"),
function (object)
{
  .oci.ResultSummary(object)
  invisible()
}
)

##
## DBIResult: Data conversion
##
setMethod("dbSetDataMappings",
signature(res = "OraResult", flds = "data.frame"),
function(res, flds, ...) .NotYetImplemented()
)

##
## DBIResult: DBI extensions
##
setGeneric("execute",
function(res, ...) standardGeneric("execute")
)

setMethod("execute",
signature(res = "OraResult"),
function(res, data = NULL, ...) .oci.execute(res, data = data)
)
