#'@title Search Guardian Tags
#'@description \code{guardian_tags} lets you directly search
#'The Guardian for their tags - metadata used to classify content.
#'
#'@param api_key A key to the Guardian API, which can be obtained
#'\href{http://open-platform.theguardian.com/access/}{here}.
#'
#'@param query Your search query. This can contain operators (\code{sausage AND mash}) or
#'phrases (\code{"sausage & mash"}); by default, searches work as an OR, looking for
#'the presence of any one individual word in the query.
#'
#'@param section the section, or sections, of \emph{The Guardian} that you want to limit the search to. Multiple
#'sections may be concatenated together using boolean operators; see \code{\link{guardian_and}} and \code{\link{guardian_and}}.
#'
#'@param reference the references to limit the search to; only tags that include those references (and meet other
#'conditions) will be returned. Also accepts boolean operators.
#'
#'@param reference_type the type of reference (such as \code{isbn}). Also accepts boolean operators.
#'
#'@param page a particular page of results to return. Useful when returning multiple sets of data with the same query;
#'you can repeat the query, incrementing the value in \code{page}.
#'
#'@param page_size the maximum number of items to return; anywhere between 1 and 50. Set to 50 by default.
#'
#'@param ... further arguments to pass to httr's \code{GET}.
#'
#'@seealso \code{\link{guardian_content}}.
#'
#'@examples
#'\dontrun{
#'# Simple example
#'results <- guardian_tags("test", "green")
#'}
#'@export
guardian_tags <- function(api_key, query, section = NULL, reference = NULL, reference_type = NULL,
                          page = NULL, page_size = 50, ...){
  
  # Base query
  path <- paste0("tags?q=", curl::curl_escape(query), "&api-key=", api_key, "&page-size=", page_size)
  
  # Filter types
  if(!is.null(section)){
    path <- paste0(path, "&section=", section)
  }
  if(!is.null(reference)){
    path <- paste0(path, "&reference=", reference)
  }
  if(!is.null(reference_type)){
    path <- paste0(path, "&reference-type=", reference_type)
  }
  
  retrieved_data <- guardian_query(path, ...)[[1]]
  
  return(retrieved_data)
}
