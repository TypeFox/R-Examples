#' Read Meta Information from a Time Series Database
#' 
#' This function reads meta information from a timeseriesdb package PostgreSQL
#' database and puts into a meta information environment. 
#' 
#' @param series character name of a time series object.
#' @param con PostgreSQL connection object
#' @param overwrite_objects logical should the entire object for a key be overwritten. Defaults to FALSE.
#' @param overwrite_elements logical should single elements inside the environment be overwritten. Defaults to TRUE.
#' @param locale character denoting the locale of the meta information that is queried.
#' defaults to 'de' for German. At the KOF Swiss Economic Institute meta information should be available
#' als in English 'en', French 'fr' and Italian 'it'. Set the locale to NULL to query unlocalized meta information. 
#' @param tbl character name of the table that contains meta information. Defaults to 'meta_data_localized'. Choose meta 'meta_data_unlocalized' when locale is set to NULL. 
#' @param meta_env environment to which the meta information should be added. Defaults to NULL. In this case an environment will be returned. If you run this function in a loop best create an empty environment before the loop or apply call and pass the environment to this function. By doing so new elements will be added to the environment. 
#' @param schema SQL schema name. Defaults to timeseries.
#' @export 
readMetaInformation <- function(series,
                                con,
                                locale = 'de',
                                tbl = 'meta_data_localized',
                                overwrite_objects = F,
                                overwrite_elements = T,
                                meta_env = NULL,
                                schema = 'timeseries'){
  
  
  tbl = paste(schema,tbl,sep='.')
  
  if(!is.null(locale)){
    sql_statement <- sprintf("SELECT (each(meta_data)).key,
                             (each(meta_data)).value
                             FROM %s WHERE ts_key = '%s' AND locale_info = '%s'",
                             tbl,series,locale)
    res <- DBI::dbGetQuery(con,sql_statement)
    res_list <- as.list(res$value)
    names(res_list) <- res$key
        
    # returns an environment of class meta_env
    addMetaInformation(series,res_list,
                       overwrite_objects =  overwrite_objects,
                       overwrite_elements = overwrite_elements,
                       meta_env = meta_env)    
  } else {
    # sanity check
    if(tbl != paste(schema,'meta_data_unlocalized',sep=".")) warning('DB table is not set to unlocalized, though locale is NULL!')
    sql_statement <- sprintf("SELECT
                             md_generated_by,
                             md_resource_last_update,
                             md_coverage_temp
                             FROM %s WHERE ts_key= '%s'",
                             tbl,series)
    
    sql_statement_hstore <- sprintf("SELECT 
                             (each(meta_data)).key,
                             (each(meta_data)).value
                             FROM %s WHERE ts_key = '%s'",tbl,series)
    res_list <- list()
    res_list$fixed <- DBI::dbGetQuery(con,sql_statement)
    res_list$flexible <- DBI::dbGetQuery(con,sql_statement_hstore)
    addMetaInformation(series,res_list,overwrite_objects =  overwrite_objects,
                       overwrite_elements = overwrite_elements,
                       meta_env = meta_env)    
  }
}

