#' Retrieve a data.frame with requested cbs tables
#' 
#' \code{get_table_list} by default a list of all tables and all columns will be retrieved.
#' You can restrict the query by supplying multiple filter statements or by specifying the
#' columns that should be returned.
#' 
#' @note \code{get_table_list} will cache results, so subsequent calls will be much faster.
#' 
#' @param ... filter statement to select rows, e.g. Language="nl"
#' @param select \code{character} columns to be returned, by default all columns
#' will be returned.
#' @param base_url optionally specify a different server. Useful for
#' third party data services implementing the same protocal.
#' @return \code{data.frame} with identifiers, titles and descriptions of tables
#' @importFrom whisker whisker.render
#' @importFrom jsonlite fromJSON
#' @export
#' @examples 
#' \dontrun{
#' 
#' # get list of english tables
#' tables_en <- get_table_list(Language="en")
#'
#' # get list of dutch tables
#' tables_nl <- get_table_list(Language="nl")
#' View(tables_nl)
#' }
get_table_list <- function(..., select=NULL, base_url = CBSOPENDATA){
  url <- whisker.render("{{BASEURL}}/{{CATALOG}}/Tables?$format=json"
                       , list( BASEURL = base_url
                             , CATALOG = CATALOG
                             )
                       )
  url <- paste0(url, get_query(..., select=select))
  
  tables <- resolve_resource(url, "Retrieving tables from")
  tables
}


## testing

# library(dplyr)
# tables <- get_table_list(Language="nl", select=c("ShortTitle","Summary"))
