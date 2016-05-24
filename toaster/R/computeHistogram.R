#' Compute histogram distribution of the column.
#' 
#' Compute histogram of the table column in Aster by mapping its value to 
#' bins based on parameters specified. When column is of numeric or
#' temporal data type it uses map-reduce histogram function over continuous 
#' values. When column is categorical (character data types) it defers to 
#' \code{\link{computeBarchart}} that uses SQL aggregate \code{COUNT(*)} with
#' \code{GROUP BY <column>}. Result is a data frame to visualize as bar charts 
#' (see creating visualizations with \code{\link{createHistogram}}).
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName Aster table name
#' @param columnName table column name to compute histogram
#' @param tableInfo pre-built summary of data to use (require when \code{test=TRUE}). See \code{\link{getTableSummary}}.
#' @param columnFrequency logical indicates to build histogram of frequencies of column
#' @param binMethod one of several methods to determine number and size of bins: \code{'manual'} indicates to use 
#'   paramters below, both \code{'Sturges'} or \code{'Scott'} will use corresponding methods of computing number
#'   of bins and width (see \url{http://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width}).
#' @param binsize size (width) of discrete intervals defining histogram (all bins are equal)
#' @param startvalue lower end (bound) of values to include in histogram
#' @param endvalue upper end (bound) of values to include in histogram
#' @param numbins number of bins to use in histogram
#' @param useIQR logical indicates use of IQR interval to compute cutoff lower and upper bounds for values to be included in 
#'   histogram: \code{[Q1 - 1.5 * IQR, Q3 + 1.5 * IQR], IQR = Q3 - Q1}
#' @param datepart field to extract from timestamp/date/time column to build histogram on
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The creteria are expressed in the form of SQL predicates (inside
#'   \code{WHERE} clause).
#' @param by for optional grouping by one or more values for faceting or alike
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \link{RODBC} 
#'   functions like \link{sqlQuery} and \link{sqlSave}).
#' @param oldStyle logical indicates if old style histogram paramters are in use (before Aster AF 5.11)
#' 
#' @export
#' @seealso \link{computeBarchart} and \link{createHistogram}
#' 
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' # Histogram of team ERA distribution: Rangers vs. Yankees in 2000s
#' h2000s = computeHistogram(channel=conn, tableName='pitching_enh', columnName='era',
#'                           binsize=0.2, startvalue=0, endvalue=10, by='teamid',
#'                           where="yearID between 2000 and 2012 and teamid in ('NYA','TEX')")
#' createHistogram(h2000s, fill='teamid', facet='teamid', 
#'                 title='TEX vs. NYY 2000-2012', xlab='ERA', ylab='count',
#'                 legendPosition='none') 
#' }
computeHistogram <- function(channel, tableName, columnName, tableInfo = NULL, 
                             columnFrequency = FALSE, binMethod = 'manual',
                             binsize = NULL, startvalue = NULL, endvalue = NULL, numbins = NULL,
                             useIQR = TRUE, datepart = NULL,
                             where = NULL, by = NULL, test = FALSE, oldStyle = FALSE) {
  
  # match argument values
  binMethod = match.arg(binMethod, c('manual', 'Sturges', 'Scott'))
  
  if (missing(tableName) || missing(columnName) || 
        is.null(tableName) || is.null(columnName)) {
    stop("Must provide table and column names.")
  }
  
  isValidConnection(channel, test)
  
  where_clause = makeWhereClause(where)
  
  if (!missing(by)) {
    byClause = paste0(ifelse(oldStyle, "BY(", "GROUP_COLUMNS("), paste0("'", by, "'", collapse=", "), ")")
    byPartition = paste(by, collapse=", ")
    bySelect = paste0(byPartition, ", ")
  }else {
    byClause = " "
    byPartition = " 1 "
    bySelect = " "
  }
  
  if (columnFrequency) {
    return (computeHistogramOfFrequencies(channel, tableName, columnName, 
                                                binsize, startvalue, endvalue, numbins,
                                                where_clause, by, byClause, byPartition, bySelect, test))
  }
  
  if (missing(tableInfo) && test) {
    stop("Must provide tableInfo when test==TRUE.")
  }
  
  if (missing(tableInfo) || !all(columnName %in% tableInfo$COLUMN_NAME)) {
    column_stats = getTableSummary(channel, tableName, include=columnName, 
                                   where=where, mock=test)
  }else {
    column_stats = tableInfo[tableInfo$COLUMN_NAME==columnName, ]
  }
    
  # check if histogram is for character column
  # if so use computeBarchart
  if (isCharacterColumn(column_stats, columnName)) {
    return (computeBarchart(channel, tableName, columnName, "COUNT(*) cnt", where=where, 
                            by=by, test=test))
  }
  
  # set number of bins default if NULL
  if (binMethod=='manual' & is.null(numbins)) {
    numbins = 30
  }
  
  # check if histogram is for date/time column
  # if so use EXTRACT function and SQL/MR
  if (isTemporalColumn(column_stats, columnName)) {
    return (computeDateHistogram(channel, tableName, columnName, 
                                 binsize, startvalue, endvalue, numbins,
                                 useIQR,
                                 where_clause, by, datepart, test))
  }
  
  # compute histogram parameters if missing
  if (binMethod=='manual') {
    if (is.null(binsize) | is.null(startvalue) | is.null(endvalue)) {
      IQR = column_stats[[1,"IQR"]]
      MIN = column_stats[[1,"minimum"]]
      MAX = column_stats[[1,"maximum"]]
      Q1 = column_stats[[1,"25%"]]
      Q3 = column_stats[[1,"75%"]]
    
      if (is.null(startvalue)) {
        if (useIQR) {
          startvalue = max(MIN, Q1-1.5*IQR)
        }else {
          startvalue = MIN
        }
      }
    
      if (is.null(endvalue)) {
        if (useIQR) {
          endvalue = min(MAX, Q3+1.5*IQR)
        }else {
          endvalue = MAX
        }
      }
      
      if (is.null(binsize)) {
        binsize = (endvalue - startvalue) / numbins
      }
    }
    
    if (startvalue >= endvalue) {
      stop("Start value should not be greater than or equal to end value. Try to run with useIQR=FALSE or check that data is not constant.")
    }
    
    histPrep = paste0(" binsize('", binsize, "')
                        startvalue('", startvalue, "')
                        endvalue('", endvalue, "') ")
    
  }else {
    histPrep = paste0(" ON hist_prep(
                          ON (SELECT ", " cast(", columnName, " as numeric) ", columnName, 
                              " FROM ", tableName, where_clause, 
                          "   )
                          VALUE_COLUMN('", columnName, "') 
                           ) as data_stat DIMENSION 
                        BIN_SELECT('", binMethod, "') ")
  }
  
  sql = paste0("SELECT * 
                  FROM hist_reduce(
                         ON hist_map(
                           ON (SELECT ", bySelect, " cast(", columnName, " as numeric) ", columnName, 
                               " FROM ", tableName, where_clause,
                           "  ) as data_input PARTITION BY ANY ",
                           histPrep,
                       "   VALUE_COLUMN('", columnName, "') ",
                           byClause,
                       "   ) 
                         partition by ", byPartition,
                      ")")
  
  
  if (test) {
    return (sql)
  }else {
    return (toaSqlQuery(channel, sql))
  }
  
}


computeHistogramOfFrequencies <- function(channel, tableName, columnName, 
                                                binsize, startvalue, endvalue, numbins,
                                                where_clause, by, byClause, byPartition, bySelect, test) {
  
  if (is.null(by)) {
     sql = paste0("SELECT * 
             FROM hist_reduce(
                    ON hist_map(
                      ON (SELECT \"", columnName, "\", count(*) cnt 
                            FROM ", tableName, where_clause,
                          " GROUP BY 1 
                        )
                      binsize('", binsize, "')
                      startvalue('", startvalue, "')
                      endvalue('", endvalue, "')
                      value_column('cnt')                                          
                    ) 
                    partition by 1 
                  )")

  }else {
     sql = paste0("SELECT * 
             FROM hist_reduce(
                    ON hist_map(                                          
                      ON (SELECT \"", by, "\", \"", columnName, "\", count(*) cnt 
                            FROM ", tableName, where_clause,
            "               GROUP BY 1, 2  
                        )
                      binsize('", binsize, "')
                      startvalue('", startvalue, "')
                      endvalue('", endvalue, "')
                      value_column('cnt')
                      by('\"", by, "\"')
                    ) 
                    partition by \"", by, "\" 
                  )")
  }
  
  if (test) {
    return (sql)
  }else {
    return (histogram = toaSqlQuery(channel,sql))
  }
  
}


computeDateHistogram <- function(channel, tableName, columnName, 
                                 binsize=NULL, startvalue=NULL, endvalue=NULL, numbins=NULL,
                                 useIQR=NULL,
                                 where_clause, by=NULL, datepart='DAY', test) {
  
  if (is.null(binsize) | is.null(startvalue) | is.null(endvalue)) {
    
    if (test) {
      stop("Must have parameters binsize, startvalue, and endvalue specified for datetime types and test==TRUE.")
    }
    
    # compute percentiles first
    percentiles = toaSqlQuery(channel,
                           paste0("SELECT * FROM approxPercentileReduce(
                                    ON (
                                      SELECT * FROM approxPercentileMap(
                                      ON  ( SELECT EXTRACT('", datepart, "' FROM \"", columnName, "\" ) dp FROM " , tableName, where_clause, " ) ",
                                  "   TARGET_COLUMN( 'dp' )
                                      ERROR( 1 )
                                      )
                                    )
                                    PARTITION BY 1
                                    PERCENTILE( 0,25,50,75,100 ))")                     
    )
    MIN = percentiles[[which(percentiles$percentile==0),"value"]]
    MAX = percentiles[[which(percentiles$percentile==100),"value"]]
    Q1 = percentiles[[which(percentiles$percentile==25),"value"]]
    Q3 = percentiles[[which(percentiles$percentile==75),"value"]]
    IQR = Q3 - Q1
    
    if (is.null(startvalue)) {
      if (useIQR) {
        startvalue = max(MIN, Q1-1.5*IQR)
      }else {
        startvalue = MIN
      }
    }
    if (is.null(endvalue)) {
      if (useIQR) {
        endvalue = min(MAX, Q3+1.5*IQR)
      }else {
        endvalue = MAX
      }
    }
    if (is.null(binsize)) { 
      binsize = (endvalue - startvalue) / numbins
    }
  }
  
  # No by clause - single histogram
  if (is.null(by)) {
       sql = paste0("SELECT * 
                       FROM hist_reduce(
                              ON hist_map(
                                ON (SELECT EXTRACT('", datepart, "' FROM \"", columnName, "\" ) dp FROM ", tableName, where_clause,
                      "  ) 
                                binsize('", binsize, "')
                                startvalue('", startvalue, "')
                                endvalue('", endvalue, "')
                                value_column('dp')
                              ) 
                              partition by 1
                            )")
                         
  # By clause - multiple histograms for each value of 'by' attribute
  }
  
  if (test) {
    return(sql)
  }else {
    return (toaSqlQuery(channel, sql))
  }
  
}