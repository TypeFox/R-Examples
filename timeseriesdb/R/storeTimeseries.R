#' Write an R Time Series to a PostgreSQL database 
#' 
#' This function writes time series object into a relational PostgreSQL database make use 
#' of PostgreSQL own 'key'=>'value' storage called hstore. The schema and database needs to 
#' created first. The parent R Package of this functions suggests a database structure
#' designed to store a larger amount of time series. This function uses INSERT INTO instead of the more convenient dbWritetable for performance reasons. DO NOT USE THIS FUNCTIONS IN LOOPS OR LAPPLY! This function can handle a set of time series on it's own and is much faster than looping over a list. Non-unique primary keys are overwritten !
#' 
#' @author Matthias Bannert, Gabriel Bucur
#' @param series character name of a time series, S3 class ts. When used with lists it is convenient to set series to names(li). Note that the series name needs to be unique in the database!
#' @param con a PostgreSQL connection object.
#' @param li list of time series. Defaults to NULL to no break legacy calls that use lookup environments.
#' @param tbl character string denoting the name of the main time series table in the PostgreSQL database.
#' @param md_unlocal character string denoting the name of the table that holds unlocalized meta information.
#' @param lookup_env environment to look in for timeseries. Defaults to .GobalEnv.
#' @param overwrite logical should existing records (same primary key) be overwritten? Defaults to TRUE.
#' @param schema SQL schema name. Defaults to timeseries. 
#' @importFrom DBI dbGetQuery
#' @export
storeTimeSeries <- function(series,
                      con,
                      li = NULL,
                      tbl="timeseries_main",
                      md_unlocal = "meta_data_unlocalized",
                      lookup_env = .GlobalEnv,
                      overwrite = T,
                      schema = "timeseries"){
  # make this function compatible with former version that used environments. 
  if(is.null(li)){
    li <- as.list.environment(lookup_env)
  }
  
  # avoid overwriting
  if(!overwrite){
    db_keys <- DBI::dbGetQuery(con,sprintf("SELECT ts_key FROM %s",tbl))$ts_key
    series <- series[!(series %in% db_keys)]
    li <- li[series]
  }
  
  # CREATE ELEMENTS AND RECORDS ---------------------------------------------
  # use the form (..record1..),(..record2..),(..recordN..)
  # to be able to store everything in one big query
  hstores <- unlist(lapply(li,createHstore))
  freqs <- sapply(li,frequency)
  values <- paste(paste0("('",
                         paste(series,
                               hstores,
                               freqs,
                               sep="','"),
                         "')"),
                  collapse = ",")
  
  
  # add schema name
  tbl <- paste(schema,tbl,sep = ".")
  md_unlocal <- paste(schema,md_unlocal,sep = ".")
  
  # CREATE META INFORMATION -------------------------------------------------
  # automatically generated meta information
  md_generated_by <- Sys.info()["user"]
  md_resource_last_update <- Sys.time()
  md_coverages <- unlist(lapply(li,function(x){
    sprintf('%s to %s',
            min(zooLikeDateConvert(x)),
            max(zooLikeDateConvert(x))
    )}
  ))
  
  # same trick as for data itself, one query
  md_values <- paste(paste0("('",
                            paste(series,
                                  md_generated_by,
                                  md_resource_last_update,
                                  md_coverages,
                                  sep="','"),
                            "')"),
                     collapse = ",")

  
  # SQL STATEMENTS ---------------------------------------------------------- 
  # data
  sql_query <- sprintf("INSERT INTO %s (ts_key,ts_data,ts_frequency) VALUES %s",
                       tbl,values) 
  # meta data
  sql_query_md <- sprintf("INSERT INTO %s (ts_key,md_generated_by,
                          md_resource_last_update,md_coverage_temp) VALUES %s",
                          md_unlocal,md_values) 
  
  # Send queries
  main_ok <- DBI::dbGetQuery(con,sql_query)
  md_ok <- DBI::dbGetQuery(con,sql_query_md)
  
  l <- length(series)
  
  if(is.null(main_ok) & is.null(md_ok)){
    paste0(l, " data and meta data records written successfully.")
  } else {
    paste("An error occured, data could not be written properly. Check the database")
  } 
}

