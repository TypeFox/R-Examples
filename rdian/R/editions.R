#'@title Search Guardian Editions
#'@description \code{guardian_editions} lets you directly search
#'The Guardian for metadata about their editions - country-specific
#'(or international) standard releases for their website's content.
#'
#'@param api_key A key to the Guardian API, which can be obtained
#'\href{http://open-platform.theguardian.com/access/}{here}.
#'
#'@param query Your search query.
#'
#'@param ... further arguments to pass to httr's \code{GET}.
#'
#'@seealso \code{\link{guardian_sections}} for retrieving sections, another category
#'of content.
#'
#'@examples
#'\dontrun{
#'# Simple example
#'uk_edition_data <- guardian_editions("test", "uk")
#'}
#'@export
guardian_editions <- function(api_key, query, ...){
  
  path <- paste0("editions?q=", curl::curl_escape(query), "&api-key=", api_key)
  retrieved_data <- guardian_query(path, ...)[[1]]
  return(retrieved_data)
}
