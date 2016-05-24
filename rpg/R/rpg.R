#' @name rpg-package
#' @aliases rpg
#' @docType package
#' @title Easy Access to Advanced PostgreSQL Features
#' @description
#' Provides functions for connecting to, reading from and writing to a PostgreSQL
#' database. Facilities for tracing the communication between R and PostgreSQL are
#' provided, as are function to retieve detailed session metadata.
#' @details
#' \tabular{ll}{
#' Package: \tab rpg\cr
#' Type: \tab Package\cr
#' Version: \tab 1.4\cr
#' Date: \tab 2015-2-21\cr
#' License: \tab GPL \cr
#' }
#' The main functions are \code{connect}, which establishes a connection,
#' \code{query}, which issues queries and \code{fetch}, which retieves results.
#' Intelligent defaults are used throughout. Functions that require a connection
#' will automatically attempt to establish a valid connection based on a previous
#' connection or from defaults. The defaults can be overriden in a variety of
#' ways.
#' @author
#' Timothy H. Keitt \cr \url{http://www.keittlab.org/} \cr \cr
#' Maintainer: Timothy H. Keitt \email{tkeitt@@gmail.com} \cr
#' @references \url{http://github.com/thk686/rpg}, \url{http://www.postgresql.org/}
#' @keywords package
#' @import Rcpp
#' @import RApiSerialize
#' @useDynLib rpg
NULL

#' PostgreSQL connection
#' 
#' Manage database connection
#' 
#' @param dbname name of the database or a valid \code{libpq} connection string
#' @param ... named optional connection parameters
#' 
#' @details If no connection parameters are supplied, the
#' connection will fallback to default parameters. Usually
#' this establishes a connection on the localhost to a database,
#' if it exists, with the same name as the user.
#'
#' Valid keywords and their defaults can be obtained by calling
#' \code{get_conn_defaults(all = TRUE)}. A valid \code{libpq}
#' connection string is composed of \code{keyword = value} pairs
#' separated by whitespace. You can either pass the entire string
#' or use named arguments. The names of the arguments will be used
#' as keywords and their values as values.
#' 
#' If a password was required but not provided, \code{connect} will
#' will open a dialog and prompt for a password. The connection is
#' then re-tried and the status returned.
#' 
#' @note Do not open a connection and then fork the R
#' process. The behavior will be unpredictable. It is perfectly
#' acceptable however to call \code{connect} within each
#' forked instance.
#' 
#' @return
#' \code{connect} returns one of:
#' \tabular{ll}{
#' \code{CONNECTION_OK} \tab Succesful connection \cr
#' \code{CONNECTION_BAD} \tab Connection failed \cr}
#' 
#' @author Timothy H. Keitt
#' 
#' @examples
#' \dontrun{
#' fetch("SHOW search_path") # default connection
#' connect("test")
#' connect(dbname = "test")
#' connect(dbname = "test", host = "localhost")
#' connect("dbname = test host = localhost")
#' disconnect()}
#' 
#' @export connect
#' @rdname connection
connect = function(dbname, ...)
{
  if ( missing(dbname) )
    values = list(...)
  else
    values = list(dbname = dbname, ...)
  if ( length(values) == 0 )
    return(connect_(character(0), character(0)))
  keywords = names(values)
  if ( is.null(keywords) || "" %in% keywords )
    stop("all arguments must be named")
  status = connect_(keywords, as.character(values))
  if ( status == "CONNECTION_BAD" &&
       get_conn_info("password.needed") &&
       interactive() )
  {
     pw = get_pw()
     keywords = c(keywords, "password")
     values = c(values, pw)
     return(connect_(keywords, as.character(values)))
  }
  return(status)
}

#' @details \code{fetch} returns the result of a query as a data frame. If
#' \code{sql} is \code{NULL} or empty, then an attempt will be made to retrieve
#' any pending resutls from previous queries. Note that query results are not
#' cleared until the next query is issued so \code{fetch} will continue to
#' return results until a new query is issued.
#' 
#' @return \code{fetch} returns a data frame or a query status object on failure.
#' @rdname query
#' @export
fetch = function(sql = "", pars = NULL)
{
  if ( is.null(sql) || nchar(sql) < 1 )
    return(fetch_dataframe())
  res = query(sql, pars)
  if ( res == "PGRES_TUPLES_OK" )
    fetch_dataframe()
  else res
}

#' @param ... list of commands to be \code{\link{paste}d} together
#' @details \code{execute} is a wrapper around \code{query}. It will raise
#' an exception if the command does not complete. Exceptions can be caught with
#' \code{\link{tryCatch}}. You cannot use a parameterized
#' query with \code{execute}. Unlike \code{query} it will \code{\link{paste}} its
#' arguments into a single string.
#' @return \code{execute} the result status string
#' @rdname query
#' @export
execute = function(...)
{
  status = query(paste(...))
  if ( status == "PGRES_BAD_RESPONSE" )
    stop("Fatal protocol error; check server")
  if ( status == "PGRES_FATAL_ERROR" )
  {
    em = query_error()
    if ( nchar(em) ) stop(em)
    stop("Fatal error")
  }
  return(status)
}

#' PostgreSQL shell
#' 
#' Run PostgreSQL's psql shell interactively
#' 
#' @param psql_opts a character string passed to the psql command
#' 
#' @details
#' The \code{psql} function repeatedly queries for input and pipes it to
#' PostgreSQL's psql command. It will terminate on \code{\\q} or empty
#' input.
#' 
#' If \code{psql_opts} is an empty string, then an attempt will be made to
#' supply suitable options based on the current connection. If there is no
#' active connection, psql will fallback to complied in defaults. If
#' \code{psql_opts} is not an empty string, then it will be passed as-is to
#' psql.
#' 
#' You can type psql's escape commands as usual. Try \code{\?}. You
#' cannot use \code{\\e} or \code{\\ef} to evoke an editor. Doing strange
#' things with \code{\!} will likely hang the R session.
#' 
#' There is no way to direclty enter a database password. If one is required,
#' you can use a \href{http://www.postgresql.org/docs/9.1/static/libpq-pgpass.html}{password file}
#'or \code{\link{set_conn_defaults}}.
#' 
#' Unfortunately it is probably impossible to enable GNU readline support
#' so for example up-arrow will recall your R commands, not the psql
#' commands entered. You can always call psql from a terminal.
#' 
#' @author Timothy H. Keitt
#'
#' @seealso \code{\link{set_default_password}}
#' 
#' @export
psql = function(psql_opts = "")
{
  psql_path = Sys.which("psql")
  if ( nchar(psql_path) == 0 ) stop("psql not found")
  psql_path = proc_psql_passwd(psql_path)
  psql_opts = proc_psql_opts(psql_opts)
  con = pipe(paste(psql_path, psql_opts))
  on.exit(close(con))
  repeat
  {
    inp = readline("psql:> ")
    if ( nchar(inp) && inp != "\\q" )
    {
      if ( inp %in% c("\\e", "\\ef") )
        cat("Cannot evoke editor\n")
      else
        writeLines(inp, con)
    }
    else break
  }
}

#' PostgreSQL database information
#' 
#' Get information about tables in a database
#' 
#' @param only.names if true, just list the table names
#' 
#' @return \code{list_tables}: a vector of table names or a data frame
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{psql}}
#' 
#' @examples
#' \dontrun{
#' system("createdb rpgtesting")
#' connect("rpgtesting")
#' begin()
#'  
#' # write data frame contents
#' data(mtcars)
#' write_table(mtcars)
#'  
#' # get some information
#' list_tables()
#' describe_table("mtcars")
#' 
#' #cleanup
#' rollback()
#' disconnect()
#' system("dropdb rpgtesting")}
#' 
#' @rdname table-info
#' @export
list_tables = function(only.names = TRUE)
{
  res = fetch("SELECT n.nspname as \"Schema\", c.relname as \"Name\",
                CASE c.relkind
                  WHEN \'r\' THEN \'table\'
                  WHEN \'v\' THEN \'view\'
                  WHEN \'i\' THEN \'index\'
                  WHEN \'s\' THEN \'special\'
                  WHEN \'f\' THEN \'foreign table\'
                END as \"Type\",
                pg_catalog.pg_get_userbyid(c.relowner) as \"Owner\"
               FROM pg_catalog.pg_class c
               LEFT JOIN pg_catalog.pg_namespace n
               ON n.oid = c.relnamespace
               WHERE c.relkind IN (\'r\',\'v\',\'f\',\'\')
               AND n.nspname <> \'pg_catalog\'
               AND n.nspname <> \'information_schema\'
               AND n.nspname !~ \'^pg_toast\'
               AND pg_catalog.pg_table_is_visible(c.oid)
               ORDER BY 1,2")
  if ( length(res) < 1 ) return(res)
  if ( inherits(res, "pq.status") ) return(res)
  if ( only.names ) return(res[[2]])
  return(res)
}

#' @param tablename the name of a PostgreSQL table
#' @param schemaname if not null, look only in this schema
#' @return \code{describe_table}: a data frame with column information
#' @rdname table-info
#' @export
describe_table = function(tablename, schemaname = NULL)
{
  sql = "SELECT
          table_schema as schema,
          table_name as table,
          column_name as column,
          ordinal_position as position,
          data_type as type,
          column_default as default
        FROM
          information_schema.columns"
  if ( is.null(schemaname) )
  {
    where = "WHERE
              table_name = $1"
    order = "ORDER BY
              table_schema, ordinal_position"
    fetch(paste(sql, where, order), tablename)
  }
  else
  {
    where = "WHERE
              table_schema = $1
             AND
              table_name = $2"
    order = "ORDER BY
              ordinal_position"
    fetch(paste(sql, where, order), c(schemaname, tablename))
  }
}

#' PostgreSQL data frame IO
#' 
#' Reads and writes table to and from database
#' 
#' @param x a data frame or something convertible to a data frame
#' @param tablename the name of the table to read from or write to
#' @param pkey a column name to use as primary key
#' @param row_names a column name to write row names
#' @param schemaname the schema name
#' @param types a list of valid PostgreSQL type names
#' @param overwrite if true, destroy existing table with the same name
#' 
#' @details
#' A table is created using the current connection. If \code{pkey} does not
#' match any column name, then a new column is created with the name given by
#' \code{pkey}. Its type will be \code{serial} and it will be set as the primary
#' key. If \code{pkey} does match an existing column name, then that column will
#' be used as the primary key. Note that \code{\link{make.unique}} will be
#' called on the column names before this matching is done. If \code{row_names}
#' is a character string, the data frame row names will be stored in a column
#' with the column name given by \code{row_names}. The \code{row_names} column
#' can also be the primary key if \code{pkey} is the same as \code{row_names}.
#' 
#' If \code{row_names} is specified when calling \code{read_table}, then the
#' resulting data frame will have row names installed from the column named
#' in \code{row_names}. Note that the column named in \code{row_names} must
#' match a column specified by \code{what}. The matching column will be removed
#' from the data frame.
#' 
#' If \code{types} is not supplied, they will be computed from the classes and
#' types of the columns of input.
#' 
#' @return \code{write_table} the final query status
#' 
#' @note The entire process is wrapped within a transcation. On failure
#' at any point, the transaction will be rolled back and the database
#' unaffected.
#' 
#' Also, \code{write_table} uses SQL \code{INSERT} statements and as such
#' will be slow for large tables. You are much better off bulk loading data
#' using the \code{COPY} command outside of \code{R}.
#' 
#' @seealso \code{\link{copy_from}}
#' 
#' @author Timothy H. Keitt
#' 
#' @examples
#' \dontrun{
#' # connect using defaults
#' system("createdb rpgtesting")
#' connect("rpgtesting")
#' begin()
#'  
#' # write data frame contents
#' data(mtcars)
#' write_table(mtcars)
#' 
#' # make "cyl" primary key (will fail unique constraint)
#' write_table(mtcars, pkey = "cyl", overwrite = TRUE)
#' 
#' # also write row names to "id"
#' write_table(mtcars, row_names = "id", overwrite = TRUE)
#' 
#' # row names as primary key
#' write_table(mtcars, row_names = "id", pkey = "id", overwrite = TRUE)
#' 
#' # default R row names and only first 3 columns
#' read_table("mtcars", what = "mpg, cyl, disp", limit = 3)
#' 
#' # row names from column "id"
#' read_table("mtcars", row_names = "id", limit = 3)
#' 
#' # get row names from primary key
#' read_table("mtcars", pkey_to_row_names = TRUE, limit = 3)
#' 
#' #cleanup
#' rollback()
#' disconnect()
#' system("dropdb rpgtesting")}
#' 
#' @rdname table-io
#' @export
write_table = function(x,
                       tablename,
                       pkey = NULL,
                       row_names = NULL,
                       schemaname = NULL,
                       types = NULL,
                       overwrite = FALSE)
{
  if ( missing(tablename) )
    tablename = deparse(substitute(x))
  x = as.data.frame(x, stringsAsFactors = FALSE)
  if ( prod(dim(x)) == 0 ) stop("Empty input")
  tableid = format_tablename(tablename, schemaname)
  x = handle_row_names(x, row_names)
  colnames = make.unique(names(x), "")
  if ( is.null(types) )
    types = sapply(x, pg_type)
  else
    types = rep(types, ncol(x))
  if ( !is.null(pkey) )
    if ( pkey %in% colnames )
    {
      i = which(colnames == pkey)
      types[i] = paste(types[i], "primary key")
      colnames = dquote_esc(colnames)
      types = paste(colnames, types)
    }
    else
    {
      types = c("serial primary key", types)
      colnames = dquote_esc(colnames)
      tabcols = c(dquote_esc(pkey), colnames)
      types = paste(tabcols, types)
    }
  else
  {
    colnames = dquote_esc(colnames)
    types = paste(colnames, types)    
  }
  types = as.csv(types)
  colnames = as.csv(colnames)
  sp = savepoint()
  on.exit(rollback(sp))
  if ( overwrite ) execute("DROP TABLE IF EXISTS", tableid)
  sql = paste("CREATE TABLE", tableid, "(", types, ")")
  status = query(sql)
  if ( status == "PGRES_COMMAND_OK" )
  {
    sqlpars = paste0("$", 1:ncol(x))
    sqlpars = as.csv(sqlpars)
    sql = paste("INSERT INTO", tableid, "(", colnames, ")")
    sql = paste(sql, "VALUES (", sqlpars, ")")
    estatus = prepare(sql)(x)
    if ( estatus == "PGRES_FATAL_ERROR" ) return(estatus)
    on.exit(commit(sp))
  }
  return(status)
}

#' @param what a vector of column names
#' @param limit only return this many rows
#' @param pkey_to_row_names if true and row_names not given, use primary key column
#' 
#' @return \code{read_table} a data frame
#' @rdname table-io
#' @export
read_table = function(tablename,
                      what = "*",
                      limit = NULL,
                      row_names = NULL,
                      schemaname = NULL,
                      pkey_to_row_names = FALSE)
{
  tablename = ifelse(is.character(tablename), tablename,
                     deparse(substitute(tablename)))
  tableid = format_tablename(tablename, schemaname)
  sql = paste("select", what, "from", tableid)
  if ( !is.null(limit) ) sql = paste(sql, "limit", limit)
  res = fetch(sql)
  if ( inherits(res, "pq.status" ) ) return(res) 
  if ( pkey_to_row_names && is.null(row_names) )
    row_names = primary_key_name(tableid)
  if ( !is.null(row_names) )
  {
    row.names(res) = res[[row_names]]
    res[[row_names]] = NULL
  }
  res
}

#' @param warn if true, \code{\link{readLines}} will issue warnings
#' @param ... passed to \code{\link{readLines}}
#' 
#' @details \code{dump_conn_trace} invokes \code{\link{readLines}} on
#' the trace file.
#' 
#' @rdname tracing
#' @export 
dump_conn_trace = function(warn = FALSE, ...)
{
  con = trace_filename()
  if ( is.null(con) || file.access(con, 4) == -1 ) return(NULL)
  out = readLines(con, warn = warn, ...)
  class(out) = "pg.trace.dump"
  return(out)
}

#' @export
print.pg.trace.dump = function(x, ...)
{
  cat(x, sep = "\n", ...)
  invisible(x)
}

#' Iterator support
#' 
#' Construct a row iterator
#' 
#' @param sql any valid query returning rows
#' @param by how many rows to return each iteration
#' @param pars optional query parameters
#' 
#' @details This function generates an interator object that can be used with
#' the \code{foreach-package}.
#' 
#' It is possible to use the
#' \code{\%dopar\%} operator as shown in the example below. You must
#' establish a connection to the database on each node and in your current
#' session because the call to \code{cursor} requires it. Note that the
#' cursor's lifetime is the current transaction block, so if anything happens
#' to the transaction or you call \code{END} or \code{ROLLBACK}, then the
#' cursor will no longer function. Apparently a named SQL cursor is visible
#' to any database session, as evidenced by the example below,
#' even though it is declared within a transaction. This is not stated
#' explicitely in the PostgreSQL documentation.
#' 
#' @note There are some reports of issues using multicore (forking) with
#' RStudio.
#'  
#' @examples
#' \dontrun{
#' # example requires foreach
#' if ( ! require(foreach, quietly = TRUE) )
#'  stop("This example requires the \'foreach\' package")
#'
#' # connect using defaults
#' system("createdb rpgtesting")
#' connect("rpgtesting")
#' begin()
#'  
#' # write data frame contents
#' data(mtcars)
#' write_table(mtcars, row_names = "id", pkey = "id", overwrite = TRUE)
#' 
#' # expand rows to columns 8 rows at a time
#' x = foreach(i = cursor("SELECT * FROM mtcars", by = 8),
#'             .combine = rbind) %do% { i$mpg }
#' print(x, digits = 2)
#'         
#' # parallel example
#' if ( require(doParallel, quietly = TRUE) )
#' {
#'  # make the cluster
#'  cl = makeCluster(2)
#'  
#'  # must connect to database on each node
#'  clusterEvalQ(cl, library(rpg))
#'  clusterEvalQ(cl, connect("rpgtesting"))
#'  clusterEvalQ(cl, begin())
#'  
#'  # setup the dopar call
#'  registerDoParallel(cl)
#'  
#'  # take column averages 4 rows at a time
#'  curs1 = cursor("SELECT * FROM mtcars", by = 4)
#'  x = foreach(i = curs1, .combine = rbind, .inorder = FALSE) %dopar%
#'  {
#'    rr = paste0(range(abbreviate(i$id)), collapse = "-")
#'    pid = get_conn_info("server.pid")
#'    j = names(i) != "id"
#'    mn = signif(apply(i[, j], 2, mean), 2)
#'    c(rows = rr, backend = pid, mn)
#'  }
#'  x = as.data.frame(x)
#'  row.names(x) = x$rows
#'  x$rows = NULL
#'  print(noquote(x))
#'  
#'  clusterEvalQ(cl, rollback())
#'  stopCluster(cl)
#' }
#' 
#' #cleanup
#' disconnect()
#' system("dropdb rpgtesting")}
#' 
#' @seealso \code{foreach}, \code{\link{rollback}}, \code{\link{query}}
#' 
#' @author Timothy H. Keitt
#' @export
cursor = function(sql, by = 1, pars = NULL)
{
  check_transaction();
  cname = unique_name();
  query(paste("DECLARE", cname, "NO SCROLL CURSOR FOR", sql), pars)
  f = function()
  {
    res = fetch(paste("FETCH", by, "FROM", cname))
    if ( inherits(res, "pq.status") ) stop(res)
    if ( length(res) < 1 ) stop("StopIteration")
    return(res)
  }
  structure(list(nextElem = f, cursor_name = cname),
            class = c('cursor', 'abstractiter', 'iter'))
}

#' Prepared queries
#'
#' Prepare and execute queries
#'
#' @param sql a valid query string
#'
#' @details \code{prepare} prepares a statement for later execution.
#' It returns a function that when called executes the prepared
#' statement. Values passed to the returned function will substituted
#' for parameters in the prepared statement. If the number of parameters
#' supplied is a
#' multiple of the number of open parameters in query prepared
#' using \code{prepare}, then the prepared query will be executed
#' repeatedly for each successive set of parameters. This repeated
#' execution loop is evaluted in C++ and so is quite fast. The
#' supplied parameter values will be coerced to a matrix of the
#' appropriate dimensions. Values passed to the function will be
#' recycled to match the number of query parameters.
#' The passed parameters will be coerced to character strings.
#'
#' @note One can use pure SQL to achieve the same result.
#'
#' It is generally a good idea to wrap \code{prepare}
#' in a transaction. If not in a transaction, you cannot rollback any updates
#' and it will be much slower as PostgreSQL initiates a transaction-per-query
#' by default.
#'
#' @return A function.
#' 
#' The function can take one argument. The values will be used
#' to fill in parameters of the prepared statement. If no argument
#' is passed, the statement will be executed without any parameters.
#'
#' @author Timothy H. Keitt
#'
#' @examples
#' \dontrun{
#' # try connecting to default database
#' system("createdb rpgtesting")
#' connect("rpgtesting")
#' begin()
#'
#' # write data frame contents
#' data(mtcars)
#' write_table(mtcars)
#'
#' # delete the rows
#' query("truncate mtcars")
#' read_table(mtcars)
#'
#' # use prepare-execute to write rows
#' pars = paste0("$", 1:11, collapse = ", ")
#' sql = paste0("INSERT INTO mtcars VALUES (", pars, ")", collapse = " ")
#' f = prepare(sql)
#' f(mtcars)
#' read_table(mtcars, limit = 5)
#'
#' # cleanup
#' rollback()
#' disconnect()
#' system("dropdb rpgtesting")}
#'
#' @rdname prepare
#' @export
prepare = function(sql)
{
  stmt = unique_statement_id()
  status = prepare_(sql, stmt)
  if ( status != "PGRES_COMMAND_OK" ) stop(query_error())
  function(x = NULL)
  {
    npars = num_prepared_params(stmt)
    if ( is.null(x) || npars == 0 )
    {
      execute("EXECUTE", stmt)
    }
    else
    {
      x = matrix(format_for_send(x), ncol = npars)
      execute_prepared_(x, stmt)
    }
  }
}

#' @param what the fields to return or all if NULL
#' @details \code{get_conn_info} returns a list containing
#' information about the current connection. For
#' readability, it will print as though it is a matrix. If
#' you want to see it as a list, try \code{unclass(get_conn_info())}.
#' 
#' If \code{length(what) == 1} then \code{get_conn_info} returns
#' a scalar
#' 
#' @return get_conn_info: a list of values
#' @export get_conn_info
#' @rdname connection-utils
get_conn_info = function(what = NULL)
{
  res = get_conn_info_()
  if ( is.null(what) ) return(res)
  if ( length(what) == 1 ) return(res[[what]])
  return(res[what])
}

#' @param ... a named list of arguments giving new defaults
#' 
#' @details \code{set_conn_defaults} sets the connection defaults by calling
#' \code{\link{Sys.setenv}} and setting the environment variable associated
#' with the connection keywords returned by \code{get_conn_defaults(all = TRUE)}.
#' These settings will only last as long as the current shell session and will
#' reset after a new login.
#' 
#' @rdname connection-utils
#' @export
set_conn_defaults = function(...)
{
  # copied from Sys.setenv
  x = list(...)
  nm = names(x)
  if (is.null(nm) || "" %in% nm) 
    stop("all arguments must be named")
  # end copy
  defs = get_conn_defaults(all = TRUE)
  for ( i in seq(along = nm) )
  {
    envvar = defs$envvar[defs$keyword == nm[i]]
    if ( length(envvar) && nchar(envvar) )
      names(x)[i] = envvar
    else
      x[i] = NULL
  }
  if ( length(x) ) do.call("Sys.setenv", x)
}

#' @param password the password
#' @details \code{set_default_password} will query for a password (if not supplied)
#' and set the \code{PGPASSWORD} environment variable accordingly. This can be used
#' with \code{\link{psql}} and \code{\link{copy_to}}.
#' @rdname connection-utils
#' @export
set_default_password = function(password = NULL)
{
  if ( is.null(password) )
    password = get_pw()
  set_conn_defaults(password = password)
  invisible()
}

#' @details \code{reset_conn_defaults} unsets all environment variables returned
#' by \code{get_conn_defaults(all = TRUE)}.
#' 
#' @rdname connection-utils
#' @export
reset_conn_defaults = function()
{
  var = get_conn_defaults(all = TRUE)$envvar
  for ( v in var )
    if ( length(v) && nchar(v) )
      Sys.unsetenv(v)
}

#' Bulk read and write
#' 
#' Read from and write to a database using COPY
#' 
#' @param what a table name or sql query string
#' @param psql_opts passed directly to the psql command line
#' 
#' @details
#' These functions use the SQL COPY command and therefore are much
#' faster than \code{\link{write_table}} and possibly
#' \code{\link{read_table}}. These functions also call PostgreSQL's
#' psql command from the command line and will fail if it is not
#' found on the search path.
#' 
#' Because these functions shell out to psql you do not need an
#' active connection. By specifying \code{psql_opts} you can connect
#' to any database without affecting the active connection. If you
#' do not specify \code{psql_opts} an attempt will be made to use
#' the active connection information. If that fails,
#' psql will use default connection settings.
#' 
#' @note These functions call \code{\link{read.csv}} and
#' \code{\link{write.csv}} and so will suffer the same bandwidth
#' limitations as those functions. I argue that is good enough.
#' There is little point in reading and writing datasets too large
#' for those functions in R. Better to bulk load using psql on
#' the command line and then use \code{\link{cursor}} to read the
#' data in small bits.
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{set_default_password}}
#' 
#' @examples
#' \dontrun{
#' # example requires hflights
#' if ( ! require(hflights, quietly = TRUE) )
#'  stop("This example requires the \'hflights\' package")
#'
#' # big dataset
#' data(hflights)
#' dim(hflights)
#' 
#' system(paste("createdb rpgtesting"))
#' 
#' opts = paste("-d rpgtesting")
#' system.time(copy_to(hflights, psql_opts = opts))
#' system.time(invisible(copy_from("hflights", psql_opts = opts)))
#' 
#' connect("rpgtesting")
#' begin()
#' 
#' ## Sloooowwwwwww
#' ## system.time(write_table(hflights))
#' system.time(invisible(read_table("hflights")))
#' 
#' rollback()
#' disconnect()
#' system(paste("dropdb rpgtesting"))}
#' 
#' @rdname copy
#' @export
copy_from = function(what, psql_opts = "")
{
  psql_path = Sys.which("psql")
  if ( nchar(psql_path) == 0 ) stop("psql not found")
  psql_opts = proc_psql_opts(psql_opts)
  if ( grepl("select", tolower(what)) ) what = paste("(", what, ")")
  sql = paste("COPY", what, "TO stdout CSV NULL \'NA\' HEADER")
  con = pipe(paste(psql_path, psql_opts, "-c", dquote_esc(sql)))
  read.csv(con, header = TRUE, as.is = TRUE)
}

#' @param x a data frame
#' @param tablename name of table to create
#' @param schemaname create table in this schema
#' @param append if false, drop and receate table
#' 
#' @rdname copy
#' @export
copy_to = function(x, tablename,
                   schemaname = NULL,
                   append = FALSE,
                   psql_opts = "")
{
  psql_path = Sys.which("psql")
  if ( nchar(psql_path) == 0 ) stop("psql not found")
  if ( missing(tablename) )
    tablename = deparse(substitute(x))
  tableid = format_tablename(tablename, schemaname)
  sql = paste("COPY", tableid, "FROM stdin CSV NULL \'NA\' HEADER")
  if ( ! append )
  {
    colnames = make.unique(names(x), "")
    types = sapply(x, pg_type)
    colspec = paste(dquote_esc(colnames), types, collapse = ", ")
    sql = paste("CREATE TABLE", tableid, "(", colspec, ");", sql)
    sql = paste("DROP TABLE IF EXISTS", tableid, ";", sql)
  }
  sql = paste("SET client_min_messages TO warning;", sql)
  psql_opts = proc_psql_opts(psql_opts)
  con = pipe(paste(psql_path, psql_opts, "-c", dquote_esc(sql)))
  write.csv(x, con, row.names = FALSE)
}

#' Transaction support
#' 
#' Start, commit or rollback transactions or savepoints
#' 
#' @details
#' These functions allow manipulation of database transaction states. If no
#' \code{savepoint} object is supplied, then an attempt is made to commit or
#' rollback the current transaction.
#' 
#' The \code{savepoint} function will initiate a transaction if one is not
#' currently active. In that case, no actual PostgreSQL savepoint will be used.
#' Rolling back the savepoint will rollback the initiated transaction. If a
#' tranaction is active, then a named savepoint will be generated. You can
#' \code{rollback} to the database state when \code{savepoint}
#' was called or \code{commit} all changes.
#' 
#' @author Timothy H. Keitt
#' 
#' @examples
#' \dontrun{
#' system("createdb rpgtesting")
#' connect("rpgtesting")
#' begin()
#' sp1 = savepoint()
#' 
#' # nest savepoints
#' sp2 = savepoint()
#' data(mtcars)
#' write_table(mtcars, "testtab", overwrite = TRUE)
#' list_tables()
#' rollback(sp2)
#' 
#' list_tables()
#' # nest savepoints
#' sp3 = savepoint()
#' sp4 = savepoint()
#' write_table(mtcars, "testtab", overwrite = TRUE)
#' commit(sp4)
#' list_tables()
#' rollback(sp3)
#' list_tables()
#' 
#' rollback(sp1)
#' rollback()
#' disconnect()
#' system("dropdb rpgtesting")}
#' 
#' @rdname transactions
#' @export
begin = function() query("BEGIN")

#' @param savepoint an object produced by \code{savepoint}
#' @rdname transactions
#' @export
commit = function(savepoint = NULL)
{
  if ( is.null(savepoint) )
    query("COMMIT")
  else
    savepoint$commit()
}

#' @rdname transactions
#' @export
rollback = function(savepoint = NULL)
{
  if ( is.null(savepoint) )
    query("ROLLBACK")
  else
    savepoint$rollback()
}

#' @return \code{savepoint}: a savepoint object
#' @rdname transactions
#' @export
savepoint = function()
{
  if ( check_transaction() )
  {
    spname = unique_name()
    execute("SAVEPOINT", spname)
    rollback = function()
      query(paste("ROLLBACK TO", spname))
    commit = function()
      query(paste("RELEASE", spname))
  }
  else
  {
    rollback = function() query("ROLLBACK")
    commit = function() query("COMMIT")
  }
  list(rollback = rollback, commit = commit)
}

#' Object storage
#' 
#' Serialize and write R objects
#' 
#' @param ... a list of objects or object names
#' @param tablename the table for storing objects
#' @param schemaname the schema to use
#' 
#' @details These functions allow one to write out any R object to a PostgreSQL
#' database and later retrieve them into any R session. The pair \code{stow} and
#' \code{retrieve} are modeled roughly as \code{\link{save}} and \code{\link{load}}.
#' 
#' The contents of \code{...} are handled specially. If a named argument is passed,
#' then the object will be stowed and retrieved with that name. A raw string will
#' not be stowed if it matches the name of any R object; the matching R object will
#' be stowed instead. A vector of strings will however be stowed as is. Object names
#' will be coerced to valid identifiers using \code{\link{make.names}}. Names must
#' be unique and not collide with any existing object name in the same table.
#' 
#' @seealso \code{\link{RApiSerialize}}
#' 
#' @author Timothy H. Keitt
#' 
#' @examples
#' \dontrun{
#' system("createdb rpgtesting")
#' connect("rpgtesting")
#' begin()
#'
#' stow("test")
#' list_stowed()
#' stow("test")
#' list_stowed()
#' stow(x = "test")
#' list_stowed()
#' x = 1
#' stow(x)
#' list_stowed()
#' stow(y = x)
#' list_stowed()
#' rm(x)
#' retrieve(".*")
#' print(test)
#' print(x)
#' print(y)
#' delete_stowed(".*")
#' data(mtcars)
#' stow(mtcars)
#' list_stowed()
#' rm(mtcars)
#' retrieve("mtcars")
#' head(mtcars)
#' 
#' rollback()
#' disconnect()
#' system("dropdb rpgtesting")}
#' 
#' @rdname stow
#' @export
stow = function(..., tablename = "rpgstow", schemaname = "rpgstow")
{
  objlist = list(...)
  objnames = names(objlist)
  objsyms = as.character(substitute(list(...)))[-1L]
  if ( is.null(objnames) ) objnames = rep("", length(objlist))
  i = which(objnames == "")
  objlist[i] = mget(objsyms[i], ifnotfound = objlist[i], envir = parent.frame())
  objnames[i] = objsyms[i]
  objnames = make.names(objnames)
  tableid = format_tablename(tablename, schemaname)
  sp = savepoint()
  on.exit(rollback(sp))
  check_stow(tablename, schemaname)
  for ( i in seq(along = objlist) )
  {
    sql = paste("INSERT INTO", tableid,
                "(objname, object) VALUES (\'", objnames[i], "\', $1)")
    status = exec_param_serialize(sql, objlist[[i]])
    if ( status == "PGRES_FATAL_ERROR" ) return(status)
  }
  on.exit(commit(sp))
  invisible()
}

#' @rdname stow
#' @export
list_stowed = function(tablename = "rpgstow", schemaname = "rpgstow")
{
  tableid = format_tablename(tablename, schemaname)
  fetch(paste("SELECT objname FROM", tableid))
}

#' @param objnames a character vector with object names or regular expressions
#' @details \code{retrieve} assigns the stowed values to the
#' stowed names in the global environment. It will overwrite any variable that
#' has the same name as a stowed object.
#' 
#' The functions \code{retrieve} and \code{delete_stowed} use regular
#' expression matching as implemented by the
#' \href{http://www.postgresql.org/docs/9.1/static/functions-matching.html}{PostgreSQL \code{~} operator}.
#' 
#' @rdname stow
#' @export
retrieve = function(objnames, tablename = "rpgstow", schemaname = "rpgstow")
{
  env = parent.frame()
  tableid = format_tablename(tablename, schemaname)
  sql = paste("SELECT * FROM", tableid, "WHERE objname ~ $1")
  for ( n in objnames )
  {
    res = fetch_stowed(sql, n)
    lapply(names(res), function(n) assign(n, res[[n]], envir = env))
  }
  invisible()
}

#' @rdname stow
#' @export
delete_stowed = function(objnames, tablename = "rpgstow", schemaname = "rpgstow")
{
  sp = savepoint(); on.exit(rollback(sp))
  tableid = format_tablename(tablename, schemaname)
  status = prepare(paste("DELETE FROM", tableid, "WHERE objname ~ $1"))
  if ( status == "PGRES_FATAL_ERROR" ) return(status)
  status = execute_prepared_(matrix(objnames, ncol = 1))
  if ( status == "PGRES_FATAL_ERROR" ) return(status)
  on.exit(commit(sp))
}

#' @param imagename a table name for stowing the session image
#' @details \code{stow_image} and \code{retrieve_image} will stow all objects in
#' the current session and retrieve them later. Note that \code{stow_image} will
#' overwrite all existing objects stowed within \code{imagename}.
#' @rdname stow
#' @export
stow_image = function(imagename = "rpgimage", schemaname = "rpgstow")
{
  sp = savepoint(); on.exit(rollback(sp))
  delete_stowed('.*', imagename, schemaname)
  args = as.list(ls(envir = globalenv(), all.names = TRUE))
  args$tablename = imagename
  args$schemaname = schemaname
  do.call("stow", args, envir = globalenv())
  on.exit(commit(sp))
}

#' @rdname stow
#' @export
retrieve_image = function(imagename = "rpgimage", schemaname = "rpgstow")
{
  args = list(objnames = '.*',
              tablename = imagename,
              schemaname = schemaname)
  do.call("retrieve", args, envir = globalenv())
}

#' @param schemaname install in this schema
#' @details \code{enable_postgis} will attempt to install the postgis
#' extension in the named schema. The default search path is altered to
#' include the new schema.
#' @rdname misc
#' @export
enable_postgis = function(schemaname = "postgis")
{
  sp = savepoint()
  on.exit(rollback(sp))
  execute("CREATE SCHEMA", schemaname)
  execute("SET search_path TO", schemaname)
  execute("CREATE EXTENSION postgis")
  execute("SET search_path TO default")
  dpath = fetch("SHOW search_path")[[1]]
  if ( ! grepl(schemaname, strsplit(dpath, ", ")) )
    execute("ALTER DATABASE", get_conn_info("dbname"),
            "SET search_path TO", dpath, ", ", schemaname)
  on.exit(commit(sp))
}
