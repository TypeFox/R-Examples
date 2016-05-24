#' Store Meta Information to the Database
#' 
#' This function stores meta information to the database for a given time series.
#' Make sure that corresponding time series had been inserted to the main table before. 
#' 
#' @param series a character name of an time series object
#' @param con a PostgreSQL connection object
#' @param tbl name of the meta information table, defaults to localized meta data: meta_data_localized. Alternatively choose meta_data_unlocalized if you are not translating meta information.
#' @param lookup_env name of the R environment in which to look for meta information objects
#' @param locale character locale fo the metainformation. Defaults to Germen 'de'. See also \code{\link{readMetaInformation}}.
#' If locale is set to NULL unlocalized meta is updated. Make sure to change tbl to 'meta_data_unlocalized'.
#' @param overwrite logical, defaults to FALSE.
#' @param quiet logical, should there be console output for every query result ? Defaults to FALSE.
#' @param schema SQL schema name, defaults to 'timeseries'.
#' @export
storeMetaInformation <- function(series,
                                 con,
                                 tbl = 'meta_data_localized',
                                 lookup_env = 'meta_data_localized',
                                 locale = 'de',
                                 overwrite = F,
                                 quiet = F,
                                 schema = 'timeseries'){
  # get an object from the meta environment
  mi <- get(series,envir = get(lookup_env))
  
  tbl <- paste(schema,tbl,sep=".")
  
  # creata a list of hstores
  
  
  # lapply a write hstore to db
  if(!is.null(locale)){
    if(overwrite){
        hstore <- createHstore(mi)
        sql_query <- sprintf("INSERT INTO %s
(ts_key,locale_info,meta_data) VALUES 
                             ('%s','%s','%s')",
                             tbl,series,locale,hstore)
    }
    } else {
      hstore <- createHstore(mi,fct = T)
      if(tbl != paste(schema,'meta_data_unlocalized',sep=".")) warning('Locale is set to NULL and tbl is not set to meta_data_unlocalized.')
      if(overwrite){
        sql_query <- sprintf("UPDATE %s SET meta_data = %s WHERE ts_key = '%s'",
                             tbl,hstore,series)
      } else{
        # use coalesce to avoid crashing when meta_data is NULL
        # concat all the hstore parts, use hstore function from PostgreSQL
        sql_query <- sprintf("UPDATE %s set meta_data = coalesce(meta_data,'')||%s WHERE ts_key = '%s'",
                             tbl,hstore,series)
      }       
      }
      
      query_return <- DBI::dbGetQuery(con,sql_query)
  
      # return proper status messages for every lang
  if(!quiet){
    if(is.null(query_return)){
      paste0(locale,' meta information successfully written.')
    } else{
      paste0(locale,' meta information fail.')
    }
  }
    
}

