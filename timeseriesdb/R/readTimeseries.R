#' Read Time Series From PostgreSQL database
#' 
#' This function reads a time series from a postgres database,
#' that uses key value pair storage (hstore), to archive time series.
#' After reading the information from the database a standard
#' R time series object of class 'ts' is built and returned. 
#' 
#' @author Matthias Bannert, Gabriel Bucur
#' @param series character vector of series names.
#' @param con a PostgreSQL connection object
#' @param tbl character string denoting the name of the view
#' containing the json records.
#' @param schema SQL schema name. Defaults to timeseries.
#' @importFrom DBI dbGetQuery
#' @importFrom RJSONIO fromJSON
#' @export
readTimeSeries <- function(series,con,tbl = "v_timeseries_json",
                     schema = "timeseries"){
  # Create a WHERE IN SQL clause, cast to json to text
  # so RS-DBI doesn't issue a warning cause it does not know
  # json
  series <- paste(paste0("'",series,"'"),collapse=",")
  sql_query <- sprintf("SELECT ts_json_records::text FROM %s.%s WHERE ts_key IN (%s)",schema,tbl,series)
  out <- DBI::dbGetQuery(con,sql_query)$ts_json_records
  
  # identify text as json and create a list of 
  # time series 
  jsn_li <- lapply(out,function(x){
    RJSONIO::fromJSON(x)
  })
  # store time series as names for the output list
  # later on
  nms <- unlist(lapply(jsn_li,"[[","ts_key"))
  
  # create time series objects
  out_li <- lapply(jsn_li,function(x){
    freq <- "[["(x,"ts_frequency")
    ts_data <- "[["(x,"ts_data")
    d <- as.Date(names(ts_data)[1])
    y <- as.numeric(format(d,"%Y"))
    p <- as.numeric(format(d,"%m"))
    if(freq == 4){
      period <- (p -1) / 3 + 1
    } else if(freq == 12){
      period <- p
    } else if(freq == 1){
      period <- NULL
    } else {
      stop("Not a standard frequency.")
    }
    
    # create the time series object but suppress the warning of creating NAs
    # when transforming text NAs to numeric NAs
    suppressWarnings(ts(as.numeric(ts_data),
                        start=c(y,period),
                        frequency = freq))
  })
  names(out_li) <- nms
  class(out_li) <- append(class(out_li),"tslist")
  return(out_li)
}
