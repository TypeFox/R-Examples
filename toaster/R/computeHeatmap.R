#' Compute 2-dimensional multi-layered matrix for heat map visualizations.
#' 
#' Compute aggregate value(s) across two category classes represented by the 
#' table columns \code{dimension1} and \code{dimension2}. Resulting data frame
#' represents 2-dimensional multi-layered matrix where each layer comprises
#' values from single aggregate. Category columns usually are of character, 
#' temporal, or discrete types. Values are aggregates computed across 
#' category columns utilizing SQL \code{GROUP BY <dimension1>, <dimension2>}. 
#' Aggregate formula may use any SQL expressions allowed with the \code{GROUP BY}
#' as defined above. Results are usually fed into \code{\link{createHeatmap}} 
#' for heat map visualizations. If defined, parameter \code{by} expands 
#' grouping columns to be used with heat maps with faceting.
#'  
#' Result represents 2-dimensional matrix with as many data layers as there were 
#' aggregates computed. Additionally more layers defined with parameter \code{by} 
#' support facets. 
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName table name
#' @param dimension1 name of the column for for heatmap x values. This value along with \code{dimension2}
#'   are x and y scales of heatmap table.
#' @param dimension2 name of the column for for heatmap y values. This value along with \code{dimension1}
#'   are x and y scales of heatmap table.
#' @param aggregates vector with SQL aggregates to compute values for heat map. Aggregate may have optional 
#'   aliases like in \code{"AVG(era) avg_era"}. Subsequently, use in \code{createHeatmap} as color 
#'   (fill), text, and threshold values for heat map cells. 
#' @param aggregateFun deprecated. Use \code{aggregates} instead.  
#' @param aggregateAlias deprecated. Use \code{aggregates} instead.
#' @param dimAsFactor logical indicates if dimensions and optional facet columns should be converted to factors.
#'   This is almost always necessary for heat maps.
#' @param withMelt logical if TRUE then uses \pkg{reshape2} \code{\link{melt}} to transform data frame
#'  with aggregate values in designated columns into a molten data frame.
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The creteria are expressed in the form of SQL predicates (inside
#'   \code{WHERE} clause).
#' @param by for optional grouping by one or more values for faceting or alike
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#' @export
#' @seealso \code{\link{createHeatmap}} 
#' @return Data frame representing 2-dimensional multi-layered matrix to use 
#'   with \code{\link{createHeatmap}}. Matrix has as many layers as there are 
#'   aggregates computed. If \code{by} defined, data frame contains multiple 
#'   matrices for each value(s) from the column(s) in \code{by} (to support facets). 
#'   When \code{withMelt TRUE} function \code{\link{melt}} applies transforming data frame
#'   and columns with aggregate values for easy casting: expands number of rows and 
#'   replaces all aggregate columns with two: \code{variable} and \code{value}.
#' 
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' hm = computeHeatmap(conn, "teams_enh", 'franchid', 'decadeid', 'avg(w) w', 
#'                     where="decadeid >= 1950")
#' hm$decadeid = factor(hm$decadeid)
#' createHeatmap(hm, 'decadeid', 'franchid', 'w')
#' 
#' # with diverging color gradient
#' hm = computeHeatmap(conn, "teams_enh", 'franchid', 'decadeid', 'avg(w-l) wl', 
#'                     where="decadeid >= 1950")
#' hm$decadeid = factor(hm$decadeid)
#' createHeatmap(hm, 'decadeid', 'franchid', 'wl', divergingColourGradient = TRUE)
#' }
computeHeatmap <- function(channel, tableName, dimension1, dimension2, 
                           aggregates = "COUNT(*) cnt",
                           aggregateFun = NULL, aggregateAlias = NULL, 
                           dimAsFactor = TRUE, withMelt = FALSE, 
                           where=NULL, by=NULL, test=FALSE) {
  
  if (missing(tableName)) {
    stop("Must have table name.")
  }
  
  if (missing(dimension1) || missing(dimension2)) {
    stop("Must have all 2 heatmap dimensions defined to compute.")
  }
  
  # validate aggregate args
  # check for deprecated parameters first
  if (!missing(aggregateFun)) {
    toa_dep("0.2.4", "\"aggregateFun\" and \"aggregateAlias\" arguments in computeHeatmap are deprecated. Use aggregates instead.")
    
    if (length(aggregateFun) != length(aggregateAlias)) 
      stop("Lengths of parameters 'aggregateFun' and 'aggregateAlias' must be the same.")
    
    aggregates = paste(aggregateFun, aggregateAlias, sep=" ")
    
  }
  
  if (is.null(aggregates) || length(aggregates) < 1)
    stop("Must have at least one aggregate defined.")
  
  isValidConnection(channel, test)
  
  where_clause = makeWhereClause(where)
  
  aggSelectList = paste(aggregates, collapse=", ")
  
  if (is.null(by)) {
     sql = paste0("SELECT ", dimension1, ", ", dimension2, ", ", aggSelectList,
                  "  FROM ", tableName, 
                  where_clause,
                  " GROUP BY 1, 2")
  }else {
     sql = paste0("SELECT ", by, ", ",  dimension1, ", ", dimension2, ", ", aggSelectList,
                  "  FROM ", tableName, 
                  where_clause,
                  " GROUP BY 1, 2, 3")  
  }
  
  if (test) {
    return (sql)
  }else {
    heatmap = toaSqlQuery(channel, sql)
  }
  
  if (dimAsFactor) {
    colNames = c(dimension1, dimension2, by)
    heatmap[, colNames] = lapply(heatmap[, colNames], FUN = function(x) {as.factor(x)})
  }
  
  if (withMelt) {
    id.vars = c(dimension1, dimension2, by)
    heatmap = melt(heatmap, id.vars=id.vars)
  }
  
  
  return (heatmap)
}