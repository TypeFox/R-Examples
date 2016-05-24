#'@title Search Guardian Sections
#'@description \code{guardian_sections} lets you directly search
#'The Guardian for metadata about sections - individual categories of
#'content, such as world or US news.
#'
#'@param api_key A key to the Guardian API, which can be obtained
#'\href{http://open-platform.theguardian.com/access/}{here}.
#'
#'@param query Your search query. This can contain operators (\code{sausage AND mash}) or
#'phrases (\code{"sausage & mash"}); by default, searches work as an OR, looking for
#'the presence of any one individual word in the query.
#'
#'@param ... further arguments to pass to httr's \code{GET}.
#'
#'@seealso \code{\link{guardian_tags}} for retrieving tags, another class
#'of content metadata.
#'
#'@examples
#'\dontrun{
#'# Simple example
#'business_sections <- guardian_sections("test", "business")
#'}
#'@export
guardian_sections <- function(api_key, query, ...){
  
  path <- paste0("sections?q=", curl::curl_escape(query), "&api-key=", api_key)
  retrieved_data <- guardian_query(path, ...)[[1]]
  return(retrieved_data)
}
