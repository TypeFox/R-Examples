.onLoad = function(libname, pkgname)
{
  connect()
}

.onUnload = function(libpath)
{
  clean_up_all()
}

#' Convert R objects to strings
#' 
#' Prepare R objects for sending to postgresql
#' 
#' @param obj any object
#' 
#' @details R objects that will be written to postgresql must be converted to
#' characters as all data is transferred to the server as text. The S3 method
#' \code{foramt_for_send} accomplishes this. It accepts any object and returns
#' a character representation.
#' 
#' You can define new conversions by supplying your own S3 override of
#' \code{format_for_send}.
#' 
#' @export
format_for_send = function(obj)
{
  UseMethod("format_for_send", obj)
}

#' @export
format_for_send.default = function(obj) as.character(obj)

#' @export
format_for_send.list = function(obj)
{
  unlist(lapply(obj, format_for_send))
}

#' @export
format_for_send.data.frame = function(obj)
{
  unlist(lapply(obj, format_for_send))
}

#' @export
format_for_send.Date = function(obj)
{
  as.character(as.POSIXlt.Date(obj))
}

pg_type = function(x)
{
  switch(class(x),
         numeric = switch(typeof(x),
                          double = "double precision",
                          integer = "integer",
                          "text"),
         integer = "integer",
         double = "double precision",
         logical = "boolean",
         Date = "date",
         "text")
}

#' @export
print.conn.info = function(x, ...)
{
  print(as.matrix(unclass(x)), ...)
}

#' @export
print.pq.error.message = function(x, ...)
{
  cat(x)
  invisible(x)
}

print.message = function(x, terminate = "\n")
{
  if ( !is.null(x) && nchar(x) > 0 )
    cat(x, terminate)
  invisible(x)
}

is_non_empty_string = function(x)
{
  if ( is.null(x) ) return(FALSE)
  if ( nchar(x) < 1 ) return(FALSE)
  return(TRUE)
}

#' @export
print.pq.status = function(x, ...)
{
  if ( getOption("verbose") ) print.message(x)
  error.message = attr(x, "error.message")
  if ( is_non_empty_string(error.message) )
    print.message(error.message)
  else
  {
    command.status = attr(x, "command.status")
    if ( is_non_empty_string(command.status) )
      print.message(command.status)
    else if ( x %in% c("PGRES_FATAL_ERROR", "BUSY", "DONE") )
      print.message(x)
  }
  invisible(x)
}

as.csv = function(...)
{
  paste0(..., collapse = ", ")
}

dquote_esc = function(...)
{
  gsub("\"+", "\"",  paste0(paste0("\"", ...), "\""))
}

format_tablename = function(tablename, schemaname = NULL)
{
  if ( is.null(schemaname) )
    dquote_esc(tablename)
  else
    paste(dquote_esc(schemaname), dquote_esc(tablename), sep = ".")
}

handle_row_names = function(a, b)
{
  if ( !is.null(b) )
  {
    a = data.frame(row.names(a), a, stringsAsFactors = FALSE)
    names(a)[1] = b
  }
  return(a)
}

primary_key_name = function(tablename)
{
  unlist(fetch("SELECT               
                  pg_attribute.attname as pkey
                FROM
                  pg_index, pg_class, pg_attribute 
                WHERE
                  pg_class.oid = $1::regclass
                AND
                  indrelid = pg_class.oid
                AND
                  pg_attribute.attrelid = pg_class.oid
                AND 
                  pg_attribute.attnum = any(pg_index.indkey)
                AND
                  indisprimary", tablename))
}

unique_name = function()
{
  dquote_esc(uuid::UUIDgenerate())
}

unique_statement_id = function()
{
  res = uuid::UUIDgenerate()
  paste0("stmt", gsub("-", "", res))
}

table_exists = function(table, schema = NULL)
{
  sql = "select count(*) > 0 from pg_catalog.pg_tables"
  if ( is.null(schema) )
    sql = paste(sql, "where tablename = $1")
  else
    sql = paste(sql, "where tablename = $1 and schemaname = $2")
  fetch(sql, strip_quotes(c(table, schema)))[[1]]
}

strip_quotes = function(x)
{
  gsub("\"", "", x)
}

proc_psql_opts = function(psql_opts = "")
{
  if ( nchar(psql_opts) == 0 && get_conn_info("status.ok") )
  {
    host = get_conn_info("host")
    dbnm = get_conn_info("dbname")
    port = get_conn_info("port")
    user = get_conn_info("user")
    if ( !is.null(host) ) psql_opts = paste(psql_opts, "-h", host)
    if ( !is.null(dbnm) ) psql_opts = paste(psql_opts, "-d", dbnm)
    if ( !is.null(port) ) psql_opts = paste(psql_opts, "-p", port)
    if ( !is.null(user) ) psql_opts = paste(psql_opts, "-U", user)
  }
  psql_opts = paste(psql_opts, "-n -q -w -1")
  return(psql_opts)
}

proc_psql_passwd = function(psql_command)
{
  if ( get_conn_info("status.ok") &&
       get_conn_info("password.supplied") )
    psql_command = paste0("PGPASSWORD=",
                          get_conn_info("password.used"),
                          " ", psql_command)
  return(psql_command)
}

run_examples = function()
{
  eval(example(ping, run.dontrun = T), globalenv())
  eval(example(connect, run.dontrun = T), globalenv())
  eval(example(query, run.dontrun = T), globalenv())
  eval(example(trace_conn, run.dontrun = T), globalenv())
  eval(example(libpq_version, run.dontrun = T), globalenv())
  eval(example(prepare, run.dontrun = T), globalenv())
  eval(example(push_conn, run.dontrun = T), globalenv())
  eval(example(async_query, run.dontrun = T), globalenv())
  eval(example(list_tables, run.dontrun = T), globalenv())
  eval(example(write_table, run.dontrun = T), globalenv())
  eval(example(cursor, run.dontrun = T), globalenv())
  eval(example(copy_to, run.dontrun = T), globalenv())
  eval(example(savepoint, run.dontrun = T), globalenv())
  invisible()
}

check_schema = function(schemaname)
{
  sql = paste("select count(*) = 0 from",
              "information_schema.schemata",
              "where schema_name = $1")
  res = fetch(sql, schemaname)
  if ( inherits(res, "pg.status") ) stop(res)
  if (  res[[1]] ) execute("create schema", dquote_esc(schemaname))
}

check_stow = function(tablename, schemaname)
{
  check_schema(schemaname)
  sql = paste("select count(*) = 0 from",
              "information_schema.tables",
              "where table_name = $1")
  res = fetch(sql, tablename)
  if ( inherits(res, "pg.status") ) stop(res)
  if ( res[[1]] )
      execute("create table", format_tablename(tablename, schemaname),
              "(objname text primary key, object bytea)")
}

get_pw = function()
{
  # This is from code floating around the internet
  tt = tcltk::tktoplevel()
  pass = tcltk::tclVar()
  tcltk::tkpack(tcltk::tklabel(tt, text = 'Password:'))
  tcltk::tkpack(tcltk::tkentry(tt, textvariable = pass, show = '*'))
  tcltk::tkpack(tcltk::tkbutton(tt, text = "Done",
                                command = function() tcltk::tkdestroy(tt)))
  tcltk::tkwait.window(tt)
  tcltk::tclvalue(pass)
}
