#' Search the Database by Keys
#' 
#' Quick handle operator to search the database by keys. All time series whose key fit 
#' the regular expression which was handed to the operator are returned in a list. 
#' 
#' @param conObj PostgreSQL Connection object.
#' @param regexp character regular expression pattern.
#' Do not set this manually, because the quick handle operator only takes two arguments.
#' use Sys.setenv to change the schema.  
#' @rdname quickHandleOps
#' @export
"%k%" <- function(conObj,regexp){
  schema <- Sys.getenv("TIMESERIESDB_SCHEMA")
  # get time series keys that suit the regexp
  sql_query <- sprintf("SELECT ts_key FROM %s.timeseries_main WHERE ts_key ~ '%s'",
                       schema,regexp)
  keys <- dbGetQuery(conObj,sql_query)$ts_key
  
  ts_list <- readTimeSeries(keys,conObj)
  return(ts_list)
}


#' Create Custom Quick Handle Operators For Unlocalized Meta Information
#' 
#' Create '%letter%' style operator to conveniently access meta data by hstore keys. 
#' This function creates a new function operator for a particular key. Name the function 
#' operator style to get the most out of it.
#' 
#' @param key character name of the key inside the hstore. 
#' @param keep_keys logical should primary time series keys be kept? Defaults to FALSE. If set to TRUE
#' dynamically created meta information handlers always use ts_key as key type no matter the key type used for the query.
#' This comes handy when storing sets of time series. 
#' @param tbl character name of the table that holds meta data. Defaults to meta_data_unlocalized. 
#' Also supports meta_data_localized
#' @param schema character database schema name. Defaults to timeseries.
#' @export
#' @rdname quickHandleOps
createMetaDataHandle <- function(key,keep_keys = FALSE, tbl = "meta_data_unlocalized",schema = "timeseries"){
  
  sql_query <- sprintf("SELECT ts_key,
                        meta_data->'%s' AS %s
                        FROM %s.%s WHERE meta_data->'%s' ~ '",
                        key,key,schema,tbl,key)
  
  sql_query <- paste0(sql_query,"%s'")
  
  # This is bit redundant but much easier to read and understand
  # in the context of creating a function dynamically and 
  # passing unevaluated expression around as strings.
  if(keep_keys) {
    fct <- sprintf("
            function(conObj,regexp){
                   sql_query <- sprintf(\"%s\",regexp)
                   
                   key_df <- dbGetQuery(conObj,sql_query)
                   
                   ts_list <- readTimeSeries(key_df$ts_key,conObj)
                   return(ts_list)
  }",sql_query,key)
  } else {
    fct <- sprintf("
            function(conObj,regexp){
                   sql_query <- sprintf(\"%s\",regexp)
                   
                   key_df <- dbGetQuery(conObj,sql_query)
                   
                   ts_list <- readTimeSeries(key_df$ts_key,conObj)
                   names(ts_list) <- key_df[match(names(ts_list), key_df$ts_key),'%s']
                   return(ts_list)
  }",sql_query,key)
  }

  
  
  eval(parse(text = fct))
}  
  




