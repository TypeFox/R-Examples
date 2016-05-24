#' Compute correlation between pairs of columns.
#' 
#' Compute global correlation between all pairs of numeric columns in table.
#' Result includes all pairwise combinations of numeric columns in the table, with 
#' optionally limiting columns to those in the parameter \code{include} or/and
#' excluding columns defined by parameter \code{except}. Limit computation 
#' on the table subset defined with \code{where}. Use \code{output='matrix'} to produce
#' results in matrix format (compatible with function \code{\link{cor}}).
#'
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName database table name
#' @param tableInfo pre-built summary of data to use (must have with \code{test=TRUE})
#' @param include a vector of column names to include. Output never contains attributes other than in the list. 
#'   When missing all columns from \code{tableInfo} included. 
#' @param except a vector of column names to exclude. Output never contains attributes from the list.
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The creteria are expressed in the form of SQL predicates (inside
#'   \code{WHERE} clause).
#' @param output Default output is a data frame of column pairs with correlation coefficient (melt format). 
#'   To return correlation matrix compatible with function \code{\link{cor}} use \code{'matrix'} .
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \link{RODBC} 
#'   functions like \link{sqlQuery} and \link{sqlSave}).
#' @return data frame with columns:
#'   \itemize{
#'     \item \emph{corr} pair of 1st and 2d columns \code{"column1:column2"}
#'     \item \emph{value} computed correlation value
#'     \item \emph{metric1} name of 1st column 
#'     \item \emph{metric2} name of 2d column
#'     \item \emph{sign} correlation value sign \code{sign(value)} (-1, 0, or 1)
#'   }
#'   Note that while number of correlations function computes is \code{choose(N, 2)}, where \code{N} is 
#'   number of table columns specified, resulting data frame contains twice as many rows by duplicating
#'   each correlation value with swaped column names (1st column to 2d and 2d to 1st positions). This 
#'   makes resulting data frame symmetrical with respect to column order in pairs and is necessary to 
#'   correctly visualize correlation matrix with \code{\link{createBubblechart}}.
#' @seealso \code{\link{createBubblechart}} and \code{\link{showData}}.
#' @export
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' cormat = computeCorrelations(channel=conn, "pitching_enh", sqlColumns(conn, "pitching_enh"), 
#'                              include = c('w','l','cg','sho','sv','ipouts','h','er','hr','bb',
#'                                          'so','baopp','era','whip','ktobb','fip'),
#'                              where = "decadeid = 2000", test=FALSE)
#' # remove duplicate correlation values (no symmetry)
#' cormat = cormat[cormat$metric1 < cormat$metric2, ]
#' createBubblechart(cormat, "metric1", "metric2", "value", label=NULL, fill="sign")
#' }
computeCorrelations <- function(channel, tableName, tableInfo, include=NULL, except=NULL, where=NULL, 
                                output=c('data.frame','matrix'), test=FALSE) {
  
  # match argument values
  output = match.arg(output, c('data.frame','matrix'))
  
  if (test & missing(tableInfo)) {
    stop("Must provide tableInfo when test==TRUE.")
  }
  
  isValidConnection(channel, test)
  
  if (missing(tableInfo)) {
    tableInfo = sqlColumns(channel, tableName)
  }
  
  columns = getNumericColumns(tableInfo, names.only=TRUE, include=include, except=except)
  
  if (is.null(columns) || length(columns) < 2) {
    stop("Must provide at least 2 numeric columns.")
  }
  
  correlations = expand.grid(columns, columns, stringsAsFactors = FALSE)
  correlations = with(correlations, correlations[Var1<Var2,])
  correlations = apply(correlations, 1, function(x) paste(x, collapse=':'))
  
  sqlmr_correlations = paste(correlations, collapse="', '")
  sql_corr_columns = paste(columns, collapse=", ")
  
  where_clause = makeWhereClause(where)
  
  
  sql = paste0("SELECT * FROM corr_reduce(
               ON corr_map(
               ON ( SELECT ", sql_corr_columns, " FROM ", tableName, where_clause, 
               " )
               columnpairs( '", sqlmr_correlations, "')
               key_name('key')
               )
               partition by key
  )")
  
  if (test) {
    return (sql)
  }else {
    rs_corrs = toaSqlQuery(channel, sql)
  }
  
  rs_corrs = cbind(rs_corrs, t(sapply(rs_corrs$corr, 
                                      FUN=function(v) unlist(strsplit(toString(v), split=":")))))
  colnames(rs_corrs)[3] = 'metric1'
  colnames(rs_corrs)[4] = 'metric2'
  
  # make data frame symmetrical
  temp = rs_corrs
  temp[,c('metric1','metric2')] = rs_corrs[,c('metric2','metric1')]
  rs_corrs = rbind(rs_corrs, temp)
  
  # produce sign column
  signs = ifelse(sign(rs_corrs$value)>0, "1", ifelse(sign(rs_corrs$value)<0, "-1", "0"))
  rs_corrs$sign = factor(signs, levels=c("-1","0","1"), ordered=TRUE)
  
  # make metric columns ordered factors with the same levels (for sorting in plots)
  rs_corrs$metric1 = factor(rs_corrs$metric1, levels=unique(rs_corrs$metric1), ordered=TRUE)
  rs_corrs$metric2 = factor(rs_corrs$metric2, levels=unique(rs_corrs$metric1), ordered=TRUE)
  
  if (output == 'matrix') {
    corrm = acast(rs_corrs[!is.nan(rs_corrs$value),], metric1~metric2, value.var='value')
    corrm = apply(corrm, c(1,2), function(x) {if(is.na(x)) 1. else x})
    return(corrm)
  }else 
    return(rs_corrs)
}
