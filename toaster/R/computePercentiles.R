#' Compute percentiles of column values.
#' 
#' Compute percentiles including boxplot quartiles across values of column 
#' \code{columnName}. Multiple sets of percentiles achieved with the
#' parameter \code{by}. Vector \code{by} may contain arbitrary number 
#' of column names: the percentiles are computed for each combination
#' of values from these columns. Remember that when using computed
#' quartiles with function \code{\link{createBoxplot}} it can utilize
#' up to 3 columns by displaying them along the x-axis and inside
#' facets.
#'   
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName Aster table name
#' @param columnName deprecated. Use vector \code{columns} instead. 
#' @param columns names of the columns to compute percentiles on
#' @param temporal logical: TRUE indicates all columns are temporal, otherwsie numerical.
#'   Temporal percentiles have 2 values: character \code{value} representing temporal
#'   percentile (date, time, timestamp or datetime) and integer \code{epoch} value 
#'   of the number of seconds since 1970-01-01 00:00:00-00 (can be negative) or for interval 
#'   values includeing \code{time}, the total number of seconds in the interval.
#' @param percentiles integer vector with percentiles to compute. Values \code{0, 25, 50, 75, 100}
#'    will always be added if omitted for numerical types, and \code{25, 50, 75, 100} for 
#'    temporal. Percentile 0 (minimum) has to be included explicitly for temporals as its 
#'    computation affects performance more than others.
#' @param by for optional grouping by one or more values for faceting or alike. 
#'   Used with \code{\link{createBoxplot}} in combination with column name for x-axis and 
#'   wrap or grid faceting.
#' @param parallel logical: enable parallel calls to Aster database. This option requires parallel 
#'   backend enabled and registered (see in examples). Parallel execution requires ODBC \code{channel} 
#'   obtained without explicit password: either with \code{\link{odbcConnect}(dsn)} or 
#'   \code{\link{odbcDriverConnect}} calls, but not with \code{\link{odbcConnect}(dsn, user, password)}.
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The creteria are expressed in the form of SQL predicates (inside
#'   \code{WHERE} clause).
#' @param nameInDataFrame name of the column in returned data frame to store table column name(s)  
#'   defined by parameter \code{columns}. \code{NULL} indicates omit this column from the data 
#'   frame (not recommended when computing percentiles for multiple columns).
#' @param stringsAsFactors logical: should columns returned as character and not excluded by \code{as.is}
#'   and not converted to anything else be converted to factors?
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \link{RODBC} 
#'   functions like \link{sqlQuery} and \link{sqlSave}).
#' @return For numeric data function returns a data frame with percentile values organized 
#'   into following columns:
#'   \itemize{
#'     \item \emph{percentile} percentile to compute (from 0 to 100): will contain all valid values 
#'       from \code{percentiles}
#'     \item \emph{value} computed percentile
#'     \item \emph{column} table column name. Override name \code{column} with parameter \code{nameInDataFrame}
#'       or omit this column all together if \code{NULL}.
#'     \item \emph{by[1], by[2], ...} in presence of parameter \code{by}, contain values of the grouping 
#'       columns for computed percentiles (optional). 
#'   }
#'   For temporal data function returns a data frame with percentile values organized 
#'   into following columns:
#'   \itemize{
#'     \item \emph{percentile} percentile to compute (from 0 to 100): will contain all valid values 
#'       from \code{percentiles}
#'     \item \emph{value} computed percentile value converted from temporal data type to its character 
#'       representation.
#'     \item \emph{epoch} corresponding to temporal percentile value epoch: for \code{date} and 
#'       \code{timestamp} values, the number of seconds since 1970-01-01 00:00:00-00 (can be negative); 
#'       for interval values include \code{time}, the total number of seconds in the interval. 
#'     \item \emph{column} table column name. Override name \code{column} with parameter \code{nameInDataFrame}
#'       or omit this column all together if \code{NULL}.
#'     \item \emph{by[1], by[2], ...} in presence of parameter \code{by}, contain values of the grouping 
#'       columns for computed percentiles (optional).
#'   }
#' 
#' @export
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' # ipouts percentiles for pitching ipouts for AL in 2000s
#' ipop = computePercentiles(conn, "pitching", "ipouts",
#'                           where = "lgid = 'AL' and yearid >= 2000")
#' 
#' # ipouts percentiles by league
#' ipopLg = computePercentiles(conn, "pitching", "ipouts", by="lgid")
#' 
#' # percentiles on temporal columns
#' playerAllDates = computePercentiles(conn, "master_enh", 
#'                     columns=c('debut','finalgame','birthdate','deathdate'),
#'                     temporal=TRUE, percentiles=c(0))
#' createBoxplot(playerAllDates, x='column', value='epoch', useIQR=TRUE, 
#'               title="Boxplots for Date columns (epoch values)", 
#'               legendPosition="none")
#' 
#' }
computePercentiles <- function(channel, tableName, columnName = NULL, columns = columnName,
                               temporal = FALSE, percentiles = c(ifelse(temporal, 5, 0),5,10,25,50,75,90,95,100), 
                               by = NULL, where = NULL, nameInDataFrame = 'column',
                               stringsAsFactors = FALSE, test = FALSE, parallel = FALSE) {
  
  if (!is.null(columnName)) {
    toa_dep("0.2.5", "\"columnName\" argument in computePercentiles is deprecated. Use \"columns\" for column names to compute percentiles on.")
  }
  
  if (missing(channel)) {
    stop("Must provide connection.")
  } 
  
  if (missing(tableName) || is.null(tableName))
    stop("Must provide table name.")
    
  if ((missing(columnName) && missing(columns)) ||
        is.null(columns) ||
        length(columns) == 0) {
    stop("Must provide at least one column name.")
  }
  
  isValidConnection(channel, test)
  
  # percentiles
  # always add 50 (median) and 25, 75 for IQR computations
  default_percentiles = c(ifelse(temporal, 25, 0), 25, 50, 75, 100)
  if (is.null(percentiles) | !is.numeric(percentiles)) {
    percentiles = unique(default_percentiles)
  }else {
    percentiles = union(percentiles, default_percentiles) 
  }
  if (any(percentiles < 0 | percentiles > 100)) {
    stop (paste("Invalid percentile value(s) passed (below 0 or above 100): ", percentiles))
  }
  percentiles = sort(percentiles)
  
  where_clause = makeWhereClause(where)
  
  percentileStr = paste(percentiles, collapse=",")
  
  if (temporal)
    return (computeTemporalPercentiles(channel, tableName, columns, 
                                       percentiles, percentileStr, by, where_clause, nameInDataFrame,
                                       stringsAsFactors, test, parallel))
  else
    return (computeNumericPercentiles(channel, tableName, columns,
                                      percentiles, percentileStr, by, where_clause, nameInDataFrame,
                                      stringsAsFactors, test, parallel))
}


computeNumericPercentiles <- function(channel, tableName, columnNames,
                                      percentiles, percentileStr, by, where_clause, nameInDataFrame,
                                      stringsAsFactors, test, parallel) {
  
  if (is.null(by)) {
    # construct column list
    partitionByList = " 1 " 
    # construct group by list by removing aliases (if any)
    groupColumnsOpt = " "
  }else {
    # construct column list
    partitionByList = paste(by, collapse=", ")
    # construct group by list by removing aliases (if any)
    groupColumnsOpt = paste(" GROUP_COLUMNS(", paste0("'", by, "'", collapse=", "), ")", sep=" ")
  }
  
  
  if (test) {
    sql = assemblePercentileSql(tableName, where_clause, columnNames[[1]], partitionByList, percentileStr, groupColumnsOpt)
    return(sql)
  }else {
    # non-functional: eliminates NOTE 'no visible binding for global variable'
    column_name = NULL
    
    if (!parallel) {
      result = foreach(column_name = columnNames, .combine='rbind', .packages=c('RODBC')) %do% {
        sql = assemblePercentileSql(tableName, where_clause, column_name, partitionByList, percentileStr, groupColumnsOpt)
        rs = toaSqlQuery(channel, sql, stringsAsFactors=stringsAsFactors)
        if (!is.null(nameInDataFrame) && nrow(rs) > 0)
          rs[, nameInDataFrame] = column_name
        rs
      }
    }else {
      result = foreach(column_name = columnNames, .combine='rbind', .packages=c('RODBC'),
                       .errorhandling='stop') %dopar% {
        sql = assemblePercentileSql(tableName, where_clause, column_name, partitionByList, percentileStr, groupColumnsOpt)
        parChan = odbcReConnect(channel)
        rs = toaSqlQuery(parChan, sql, stringsAsFactors=stringsAsFactors, closeOnError=TRUE)
        close(parChan)
        if (!is.null(nameInDataFrame) && nrow(rs) > 0)
          rs[, nameInDataFrame] = column_name
        rs
      }
    }
    
    return(result)
  }
}

assemblePercentileSql <- function(tableName, where_clause, name, partitionByList, percentileStr, groupColumnsOpt) {
  
  sql = paste0("SELECT * FROM approxPercentileReduce(
                                  ON (
                                    SELECT * FROM approxPercentileMap(
                                      ON  ( SELECT * FROM " , tableName, where_clause, " ) ",
               "                      TARGET_COLUMN( '", name, "' )
                                      ERROR( 1 ) ",
               groupColumnsOpt,
               "                     )
                                     )
                                  PARTITION BY ", partitionByList, 
               "                  PERCENTILE( ", percentileStr, " )",
               groupColumnsOpt,
               "                  )")
}


computeTemporalPercentiles <- function(channel, tableName, columnNames,
                                       percentiles, percentileStr, by, where_clause, nameInDataFrame,
                                       stringsAsFactors, test, parallel) {
  
  if (is.null(by)) {
    # construct group select list
    selectByList = " "
    # construct partition by list
    partitionByList = " " 
    # construct group by list 
    groupByList = "1"
    groupBy0PercentileList = ""
  }else {
    # construct group select list
    selectByList = paste0(paste(by, collapse=", "), ", ")
    # construct partition by list
    partitionByList = paste("PARTITION BY", paste(by, collapse=", "))
    # construct group by list by removing aliases (if any)
    groupByList = paste(as.character(1:(length(by)+1)), collapse=", ")
    groupBy0PercentileList = paste("GROUP BY", paste(as.character(1:length(by)), collapse=", "))
  }
  
  if (test) {
    sql = assembleTemporalPercentileSql(tableName, where_clause, percentiles, percentileStr, 
                                        selectByList, partitionByList, groupByList, 
                                        groupBy0PercentileList)
    return(gsub('(%%%column__name%%%)', columnNames[[1]], sql, fixed=TRUE))
  }else {
    # non-functional: eliminates NOTE 'no visible binding for global variable'
    column_name = NULL
    
    if (!parallel) {
      result = foreach(column_name=columnNames, .combine='rbind', .packages=c('RODBC')) %do% {
        sql = assembleTemporalPercentileSql(tableName, where_clause, percentiles, percentileStr,
                                            selectByList, partitionByList, groupByList,
                                            groupBy0PercentileList)
        rs = toaSqlQuery(channel, gsub('(%%%column__name%%%)', column_name, sql, fixed=TRUE), stringsAsFactors=stringsAsFactors)
        if (!is.null(nameInDataFrame) && nrow(rs) > 0)
          rs[, nameInDataFrame] = column_name
        rs
      }
    }else {
      result = foreach(column_name = columnNames, .combine='rbind', .packages=c('RODBC'),
                       .errorhandling='stop') %dopar% {
        sql = assembleTemporalPercentileSql(tableName, where_clause, percentiles, percentileStr,
                                            selectByList, partitionByList, groupByList,
                                            groupBy0PercentileList)
        parChan = odbcReConnect(channel)
        rs = toaSqlQuery(parChan, gsub('(%%%column__name%%%)', column_name, sql, fixed=TRUE), stringsAsFactors=stringsAsFactors)
        close(parChan)
        if (!is.null(nameInDataFrame) && nrow(rs) > 0)
          rs[, nameInDataFrame] = column_name
        rs
      }
    }
    
    return (result)
  }
  
}

assembleTemporalPercentileSql <- function(tableName, where_clause, percentiles, percentileStr, 
                                          selectByList, partitionByList, groupByList, 
                                          groupBy0PercentileList) {
  
  where_clause = ifelse(grepl('[^ ]',where_clause), paste0(where_clause, " AND (%%%column__name%%%) IS NOT NULL"),
                        " WHERE (%%%column__name%%%) IS NOT NULL")
  
  sql = paste0("SELECT ", selectByList, " percentile, 
                       MAX((%%%column__name%%%))::varchar value, 
                       MAX(EXTRACT('EPOCH' FROM (%%%column__name%%%))) epoch  
              FROM (SELECT ", selectByList, " (%%%column__name%%%), 
                              ntile(100) OVER (", partitionByList, " ORDER BY (%%%column__name%%%)) percentile
                      FROM ", tableName, where_clause, ") t
             WHERE percentile IN ( ", percentileStr, " ) 
             GROUP BY ", groupByList, 
               # performance enhancement: include SELECT below only when 0 percentile (min) requested
               ifelse(0 %in% percentiles, paste0(
          " UNION  
            SELECT ", selectByList, " 0, MIN((%%%column__name%%%))::varchar, MIN(EXTRACT('EPOCH' FROM (%%%column__name%%%))) epoch 
              FROM ", tableName, where_clause, " ", groupBy0PercentileList),
                 " "),
               "  ORDER BY ", groupByList)
  
}