#' List Aster numeric data types.
#' 
#' @return character vector with names of Aster numeric data types
#' @export
#' @seealso \code{\link{getCharacterTypes}}, \code{\link{getTemporalTypes}}, \code{\link{getTableSummary}}
#' @examples
#' getNumericTypes()
#' 
getNumericTypes <- function () {
  return( c('integer',
            'numeric',
            'bigint',
            'smallint',
            'real',
            'double precision',
            'serial',
            'bigserial',
            'float',
            'decimal')
  )
}

#' List Aster character data types.
#'
#' @return character vector with names of Aster character data types
#' @export
#' @seealso \code{\link{getNumericTypes}}, \code{\link{getTemporalTypes}}, \code{\link{getTableSummary}}
#' @examples 
#' getCharacterTypes()
#' 
#' 
getCharacterTypes <- function() {
  return(c('varchar',
           'char',
           'character')
  )
}

#' List Aster temporal data types.
#' 
#' @return character vector with names of Aster temporal data types
#' @export
#' @seealso \code{\link{getCharacterTypes}}, \code{\link{getCharacterTypes}}, \code{\link{getTableSummary}}
#' @examples
#' getTemporalTypes()
#' 
getTemporalTypes <- function() {
  return(c('date', 
           'timestamp without time zone', 
           'timestamp with time zone',
           'time without time zone',
           'time with time zone')
  )
}

getTypes <- function(types) {
  
  result = character(0)
  
  if ('numeric' %in% types) {
    result = union(result, getNumericTypes())
  }
  
  if ('character' %in% types) {
    result = union(result, getCharacterTypes())
  }
  
  if ('temporal' %in% types) {
    result = union(result, getTemporalTypes())
  }
  
  return(result)
}

#' Filter numeric columns.
#'
#' Select numeric columns (names or rows) from table info data frame.
#' 
#' @param tableInfo data frame obtained by calling \code{\link{getTableSummary}}.
#' @param names.only logical: if TRUE returns column names only, otherwise full rows of \code{tableInfo}.
#' @param include a vector of column names to include. Output is restricted to this list.
#' @param except a vector of column names to exclude. Output never contains names from this list.
#' 
#' @seealso \code{\link{getCharacterColumns}}, \code{\link{getTemporalColumns}}, \code{\link{getTableSummary}}
#' @export
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' pitchingInfo = getTableSummary(channel=conn, 'pitching_enh')
#' getNumericColumns(pitchingInfo)
#' num_cols_df = getNumericColumns(pitchingInfo, names.only=FALSE)
#' }
getNumericColumns <- function (tableInfo, names.only=TRUE, include=NULL, except=NULL) {
  
  numeric_types = getNumericTypes()
  
  return(getColumns(tableInfo, numeric_types, names.only, include, except))
}



#' Filter character columns.
#' 
#' Selects character columns (names or rows) from table info data frame.
#' 
#' @param tableInfo data frame obtained by calling \code{\link{getTableSummary}}.
#' @param include a vector of column names to include. Output is restricted to this list.
#' @param except a vector of column names to exclude. Output never contains names from this list.
#' @param names.only logical: if TRUE returns column names only, otherwise full rows of \code{tableInfo}.
#' @seealso \code{\link{getNumericColumns}}, \code{\link{getTemporalColumns}}, \code{\link{getTableSummary}}
#' @export
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' pitchingInfo = getTableSummary(channel=conn, 'pitching_enh')
#' getCharacterColumns(pitchingInfo)
#' char_cols_df = getCharacterColumns(pitchingInfo, names.only=FALSE)
#' }
getCharacterColumns <- function (tableInfo, names.only=TRUE, include=NULL, except=NULL) {
  
  char_types = getCharacterTypes()
  
  return(getColumns(tableInfo, char_types, names.only, include, except))
}


#' Filter Date and Time Table Columns.
#' 
#' Selects date and time columns (names or rows) from table info data frame.
#' 
#' @param tableInfo data frame obtained by calling \code{\link{getTableSummary}}.
#' @param include a vector of column names to include. Output is restricted to this list.
#' @param except a vector of column names to exclude. Output never contains names from this list.
#' @param names.only logical: if TRUE returns column names only, otherwise full rows of \code{tableInfo}.
#' @seealso \code{\link{getCharacterColumns}}, \code{\link{getNumericColumns}}, \code{\link{getTableSummary}}
#' @export
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' masterInfo = getTableSummary(channel=conn, 'master')
#' getTemporalColumns(masterInfo)
#' date_cols_df = getTemporalColumns(masterInfo, names.only=FALSE)
#' }
getTemporalColumns <- function (tableInfo, names.only=TRUE, include=NULL, except=NULL) {
  datetime_types = getTemporalTypes()
  
  return(getColumns(tableInfo, datetime_types, names.only, include, except))
}


#' Filter columns by pattern.
#' 
#' Selects columns with names matching regular expression pattern.
#' 
#' @param pattern character string containing a \link{regular expression} to be matched in the given table info.
#' @param channel connection object as returned by \code{\link{odbcConnect}}. Only used in combination with \code{tableName}.
#' @param tableName Aster table name to use. If missing then \code{tableInfo} will be used instead.
#' @param tableInfo data frame obtained by calling \code{\link{getTableSummary}} or \code{\link{sqlColumns}}.
#' @param names.only logical: if TRUE returns column names only, otherwise full rows of \code{tableInfo}.
#' @param ignore.case if TRUE case is ignored during matching, otherwise matching is case sensitive.
#' @param invert logical. if TRUE return columns that do not match.
#' @seealso \code{\link{grep}}, \code{\link{getTableSummary}}
#' @export
getMatchingColumns <- function (pattern, channel, tableName, tableInfo, names.only = TRUE, 
                                ignore.case = TRUE, invert = FALSE) {
  
  if (!missing(tableName)) {
    tableInfo = sqlColumns(channel, tableName)
  }
  idx = grep(pattern, tableInfo$COLUMN_NAME, ignore.case=ignore.case, value=FALSE, invert=invert)
  
  if (names.only) 
    return(tableInfo[idx, "COLUMN_NAME"])
  else
    return(tableInfo[idx, ])
}


isCharacterColumn <- function (tableInfo, columnName) {
  is_column_char = getCharacterColumns(tableInfo, names.only=TRUE, include=columnName)
  return (ifelse(length(is_column_char) == 1, TRUE, FALSE))
}


isNumericColumn <- function (tableInfo, columnName) {
  is_column_numeric = getNumericColumns(tableInfo, names.only=TRUE, include=columnName)
  return (ifelse(length(is_column_numeric) == 1, TRUE, FALSE))
}


isTemporalColumn <- function (tableInfo, columnName) {
  is_column_datetime = getTemporalColumns(tableInfo, names.only=TRUE, include=columnName)
  return (ifelse(length(is_column_datetime) == 1, TRUE, FALSE))
}


includeExcludeColumns <- function (tableInfo, include, except) {
  result = tableInfo
  
  if(!is.null(include))
    result = result[result$COLUMN_NAME %in% include,]
  
  if(!is.null(except)) 
    result = result[!result$COLUMN_NAME %in% except,]
  
  return(result)
}

 
getColumns <- function (tableInfo, types, names.only, include, except) {
  result = tableInfo[tableInfo$TYPE_NAME %in% types,]
  
  result = includeExcludeColumns(result, include, except)
  
  if (names.only) 
    return(result[,"COLUMN_NAME"])
  else
    return(result)
}


makeSqlColumnList <- function(columns) {
  
  paste(columns, collapse=", ")
}


makeSqlMrColumnList <- function(columns) {
  
  paste0("'", paste(columns, collapse="', '"), "'")
}


makeSqlValueList <- function(values) {
  
  if(is.numeric(values))
    paste0(values, collapse = ", ")
  else if(is.character(values))
    paste0("'", paste(values, collapse = "', '"), "'")
  else
    stop("Values must be either numeric or character only.")
  
}


makeSqlAggregateColumnList <- function(columns, sqlAggFun, includeFunInAlias=TRUE, cast="") {
  
  if (includeFunInAlias)
    paste0(toupper(sqlAggFun), "(", columns, ")", cast, " ", tolower(sqlAggFun), '_', columns, collapse = ", ")
    # paste0(sqlAggFun, "(", columns, ")", cast, " ", paste(sqlAggFun, columns, sep='_'), collapse = ", ")
  else 
    paste0(toupper(sqlAggFun), "(", columns, ")", cast, " ", columns, collapse = ", ")
}


makeWhereClause <- function (where) {
  
  if(is.null(where))
    where_clause = " "
  else
    where_clause = paste(" WHERE", where, " ")
  
  return(where_clause)
}


makeOrderByClause <- function (order) {
  if (is.null(order))
    orderby_clause = " "
  else
    orderby_clause = paste(" ORDER BY", paste(order, collapse=", "))
  
  return (orderby_clause)
}


makeLimitClause <- function (top) {
  if (is.null(top)) 
    limit_clause = " "
  else
    limit_clause = paste(" LIMIT", top)
  
  return (limit_clause)
}

normalizeTableName <- function (name) {
  
  tolower(name)
}

#' Make Aster temporary table name.
#'
#' @param prefix Table name will always start with toa_temp_ followed by prefix (if exists).
#' @param n non-negative integer giving number of random characters to include in the name.
#' @param schema Aster database schema table should belong to. 
#' 
#' Table name generated will always begin with 'toa_temp_' followed by prefix (if not NULL) 
#' and n random alpha-numeric characters. Beware that total length can not exceed than 63 (Aster 
#' limit on table name length).
#' 
#' @return character string suitable for Aster temporary table name
#' @export
#' @seealso \code{\link{getTableSummary}}
#' @examples 
#' tempTableName = makeTempTableName("centroids", 20)
#' 
#'
makeTempTableName <- function(prefix=NULL, n=20, schema=NULL) {
  
  if(!is.null(prefix) && !grepl("^[a-z0-9]+$", prefix, ignore.case=TRUE))
    stop("Prefix may contain alpha-numeric characters only")
  
  prefix = paste0("toa_temp_", prefix, ifelse(is.null(prefix), "", "_"))
  if (nchar(prefix) + n > 63)
    stop("Too long prefix: 63 characters is Aster limit on table name length")
  
  if(!is.null(schema) && !grepl("^[a-z0-9]+$", schema, ignore.case=TRUE))
    stop("Schema may contain alpha-numeric characters only")
  
  schema = ifelse(is.null(schema), "", paste0(schema,"."))
  return(paste0(schema, prefix, paste0(sample(c(letters,0:9), n-length(prefix), replace=TRUE), collapse="")))
}


#' Make SQL FROM clause 
#' 
#' @param name table or view name or a SQL query.
#' @param flag logical indicates if a table or a query is visible.
#'   Special value \code{NA} indicates that \code{name} is
#'   a SQL query.
#' @param alias query alias to use. Ignored if \code{name} is 
#'   a table or a view.
#'  
makeFromClause <- function(name, flag, alias = 't') {
  
  if (is.null(name)) 
    stop("Table name or query is NULL.")
  
  if (is.null(flag)) 
    stop("")
  
  if(is.na(flag))
    paste0("(", name, ")", ifelse(is.null(alias), "", paste0(" ", alias)))
  else if(flag)
    name
  else
    stop(paste0("Table ", name, " not found"))
  
}



#' Determine window function to use
#' 
#' @param rankFunction one of rank function codes to map to one of SQL window
#'   functions for ranking.
#'
getWindowFunction <- function(rankFunction) {
  windowFunction = switch(tolower(rankFunction),
                          rank="RANK()",
                          row="ROW_NUMBER()",
                          rownumber="ROW_NUMBER()",
                          denserank="DENSE_RANK()",
                          percentrank = "PERCENT_RANK()"
  )
  
  return(windowFunction)
}