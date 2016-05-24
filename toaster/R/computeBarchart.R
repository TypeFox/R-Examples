#' Compute one or more aggregates across single class.
#' 
#' Compute aggregates across category class represented by the table 
#' column. Values are one or more SQL aggregates that are valid expressions
#' with \code{GROUP BY <class column>}. Class column usually is of character 
#' or other discrete type. Typical example is computing a bar chart for the 
#' column using SQL \code{COUNT(*) ... GROUP BY} - hence the name of the 
#' function. Result is a data frame to visualize as bar charts or heatmaps
#' (see creating visualizations with \code{\link{createHistogram}} and
#' \code{\link{createHeatmap}}).
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName table name
#' @param category column name or expression associated with categories. Name may be 
#'   valid SQL expression and can contain otional alias (e.g. \code{"UPPER(car_make) make"})
#' @param aggregates SQL aggregates to compute. Each aggregate corresponds to category value. 
#'   Aggregates may have optional aliases like in \code{"AVG(era) era"}
#' @param by for optional grouping by one or more columns for faceting or alike (effectively these elements
#'   will be part of \code{GROUP BY ...}) 
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The creteria are expressed in the form of SQL predicates (inside
#'   \code{WHERE} clause).
#' @param orderBy list of column names, aliases, references or their combinations to use in SQL \code{ORDER BY} 
#'   clause. Use in combination with \code{top} below to compute only limited number of results in certain order.
#' @param top if specified indicates number of bars to include in bar plot. In combination with \code{orderBy} 
#'   it works as computing first \code{top} results in certain order.
#' @param withMelt logical if TRUE then uses \pkg{reshape2} \code{\link{melt}} to transform result data frame
#'  aggregate values into a molten data frame
#' @param stringsAsFactors logical: should columns returned as character and not excluded by \code{as.is}
#'   and not converted to anything else be converted to factors?
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \link{RODBC} 
#'   functions like \link{sqlQuery} and \link{sqlSave}).
#' @return Data frame to use for bar chart plots with \code{\link{createHistogram}}.
#' @export
#' @seealso \code{\link{computeHistogram}}, \code{\link{createHistogram}}
#' 
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' # Compute average team season era, walks, and hits for each decade starting with 1980
#' computeBarchart(channel=conn, "teams_enh", "teamid team", 
#'                 aggregates=c("avg(era) era", "avg(bb) bb", "avg(h) h"), 
#'                 where="yearid >=1980", by=c("decadeid"))
#'                 
#' # multipe aggregates in the same bar chart (with melt)
#' bc = computeBarchart(channel=conn, tableName="pitching_enh", category="teamid",
#'                     aggregates=c("AVG(era) era", "AVG(whip) whip"), withMelt=TRUE,
#'                     where="yearid >= 2000 and lgid='AL'")
#'         
#' # adding facets by decadeid          
#' bc = computeBarchart(channel=conn, tablelName="pitching_enh", category="teamid", 
#'                      aggregates=c("AVG(era) era", "AVG(whip) whip", "AVG(ktobb) ktobb"),
#'                      where="yearid >= 1990 and lgid='AL'", by="decadeid", withMelt=TRUE)
#' }
computeBarchart <- function(channel, tableName, category,
                            aggregates = "COUNT(*) cnt", 
                            where = NULL, orderBy = NULL, top = NULL, by = NULL, 
                            withMelt = FALSE,
                            stringsAsFactors = FALSE, test = FALSE) {
  
  if (missing(tableName)) {
    stop("Must have table name.")
  }
  
  if (missing(category) || length(category) != 1) {
    stop("Bar chart must have exactly one category.")
  }
  
  if (is.null(aggregates) || length(aggregates) < 1) {
    stop("Must have at least one aggregate defined.")
  }
  
  isValidConnection(channel, test)
  
  where_clause = makeWhereClause(where)
  
  orderby_clause = makeOrderByClause(orderBy)
  
  limit_clause = makeLimitClause(top)
  
  categoryExpr = sub(category, pattern = " [a-zA-Z0-9_]*$", replacement = "")
  
  if (missing(by)) {
    # construct column list
    columnList = paste(category, paste(aggregates, collapse=", "), sep=", ")
    # construct group by list by removing aliases (if any)
    groupByList = paste(categoryExpr, collapse=", ")
  }else {
    # construct column list
    columnList = paste(category, paste(by, collapse=", "), paste(aggregates, collapse=", "), sep=", ")
    # construct group by list by removing aliases (if any)
    groupByList = paste(categoryExpr, paste(by, collapse=", "), sep=", ")
  }
  
  sql = paste0("SELECT ", columnList, " FROM ", tableName,  
               where_clause,
               " GROUP BY ", groupByList,
               orderby_clause, 
               limit_clause)
  
  if (test) {
    return(sql)
  }else {
    df = toaSqlQuery(channel, sql, stringsAsFactors=stringsAsFactors)
    
    if (withMelt) {
      df = melt(df, id.vars=c(category, by))
    }
    
    return(df)
  }

}