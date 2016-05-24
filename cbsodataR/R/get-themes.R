#' Get list of all cbs thematic entries.
#' 
#' Returns a list of all cbs themes. 
#' @param ... Use this to add a filter to the query e.g. \code{get_themes(ID=10)}.  
#' @param select \code{character} vector with names of wanted properties. default is all
#' @param base_url optionally specify a different server. Useful for
#' third party data services implementing the same protocal.
#' @return A \code{data.frame} with various properties of SN/CBS themes.
#' 
#' The filter is specified with \code{<column_name> = <values>} in which \code{<values>} is a character vector.
#' Rows with values that are not part of the character vector are not returned.
#' @export
#' @examples 
#' \dontrun{
#' # get list of all themes
#' get_themes()
#' 
#' # get list of all dutch themes from the Catalog "CBS"
#' get_themes(Language="nl", Catalog="CBS")
#' }
#' @importFrom whisker whisker.render
get_themes <- function(..., select=NULL, base_url = CBSOPENDATA){
  url <- whisker.render("{{BASEURL}}/{{CATALOG}}/Themes?$format=json"
                       , list( BASEURL = base_url
                             , CATALOG = CATALOG
                             )
                       )
  url <- paste0(url, get_query(..., select=select))  
  themes <- resolve_resource(url, "Retrieving themes from ")
  themes
}



## testing

# library(dplyr)
#themes <- get_themes()
