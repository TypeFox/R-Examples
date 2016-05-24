#' Delete Time Series from the database
#' 
#' This function deletes time series AND their metainformation from the 
#' database. All meta information in all series will be deleted. 
#' To only edit the original time series use \code{\link{storeTimeSeries}}
#' to overwrite existing series. 
#' 
#' @param series character name of the timeseries
#' @param con a PostgreSQL connection object
#' @param tbl_main character name of the table that contains the 
#' main time series catalog. Defaults to 'timeseries_main'.
#' @param schema SQL schema name. Defaults to 'timeseries'.
#' @export
deleteTimeSeries <- function(series,con,
                             tbl_main = 'timeseries_main',
                             schema = 'timeseries'){
  sql_statement <- sprintf("DELETE FROM %s WHERE ts_key ='%s' CASCADE",
                           tbl_main, series)
  if(is.null(DBI::dbGetQuery(con,sql_statement))) sprintf('Time series %s deleted.',series)
  
}