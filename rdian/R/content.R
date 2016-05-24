# Collapsing the values of content calls.
collapse_content <- function(api_response){
  
  response_names <- names(api_response$results[[1]])
  results <- do.call("rbind", lapply(api_response$results, function(x){
    unlisted_data <- unname(unlist(x))
    return(as.data.frame(matrix(unlisted_data, nrow = 1, ncol = length(unlisted_data)), stringsAsFactors = FALSE))
  }))
  names(results) <- response_names
  if("webPublicationDate" %in% response_names){
    results$webPublicationDate <- as.POSIXct(gsub(x = results$webPublicationDate, pattern = "(T|%Z)",
                                                  replacement = ""), tz = "UTC")
  }
  return(results)
}


#'@title Search Guardian Content
#'@description \code{guardian_content} lets you directly search
#'The Guardian for content, including support for conditional
#'searches, date-range based filtering, and section or tag
#'based filtering.
#'
#'@param api_key A key to the Guardian API, which can be obtained
#'\href{http://open-platform.theguardian.com/access/}{here}.
#'
#'@param query Your search query. This can contain operators (\code{sausage AND mash}) or
#'phrases (\code{"sausage & mash"}); by default, searches work as an OR, looking for
#'the presence of any one individual word in the query.
#'
#'@param from the date to limit by. If set, the search will only cover data written after this date.
#'Optional (NULL) by default. Can be set in conjunction with \code{to}.
#'
#'@param to Another date to limit by. If set, the search will only cover data written before this date.
#'Optional (NULL) by default. Can be set in conjunction with \code{from}.
#'
#'@param section the section, or sections, of \emph{The Guardian} that you want to limit the search to. Multiple
#'sections may be concatenated together using boolean operators; see \code{\link{guardian_and}} and \code{\link{guardian_or}}.
#'
#'@param reference the references to limit the search to; only articles that include those references (and meet other
#'conditions) will be returned. Also accepts boolean operators.
#'
#'@param reference_type the type of reference (such as \code{isbn}). Also accepts boolean operators.
#'
#'@param tags the tags to limit the search to; only articles that contain these tags will be returned. Also accepts
#'boolean operators.
#'
#'@param rights limit a search to content with these rights (syndicatable content, for example). Options are
#'\code{syndicatable} and \code{subscription-database}. Does not accept boolean operators.
#'
#'@param ids limit a search to content with these IDs, an ID being the standard URL fragment after \code{guardian.com} (
#'or \code{co.uk}, or...). Does not accept boolean operators.
#'
#'@param production_office the production office(s) to limit to, such as "aus". Accepts boolean operators.
#'
#'@param page a particular page of results to return. Useful when returning multiple sets of data with the same query;
#'you can repeat the query, incrementing the value in \code{page}.
#'
#'@param page_size the maximum number of items to return; anywhere between 1 and 50. Set to 50 by default.
#'
#'@param fields additional fields to include in the returned data (when available). Possible values
#'(which should appear in a vector) appear \href{http://open-platform.theguardian.com/documentation/search}{here}.
#'
#'@param collapse whether to collapse the actual content into a data.frame. Set to TRUE by default.
#'
#'@param ... further arguments to pass to httr's \code{GET}.
#'
#'@seealso \code{\link{guardian_tags}} to retrieve metadata about tags.
#'
#'@examples
#'
#'# Simple example
#'\dontrun{
#'result <- guardian_content("test", "debate AND economy")
#'}
#'@export
#'@importFrom curl curl_escape
guardian_content <- function(api_key, query, from = NULL, to = NULL, section = NULL,
                             reference = NULL, reference_type = NULL, tags = NULL,
                             rights = NULL, ids = NULL, production_office = NULL, page = NULL,
                             page_size = 50, fields = NULL, collapse = TRUE, ...){
  
  # Construct basic path
  path <- paste0("search?q=", curl::curl_escape(query), "&api-key=", api_key, "&page-size=", page_size)
  
  # Check dates
  if(!is.null(from)){
    path <- paste0(path, "&from-date=", date_convert(from))
  }
  if(!is.null(to)){
    path <- paste0(path, "&to-date=", date_convert(to))
  }
  
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
  if(!is.null(tags)){
    path <- paste0(path, "&tag=", tags)
  }
  if(!is.null(rights)){
    path <- paste0(path, "&rights=", rights)
  }
  if(!is.null(ids)){
    path <- paste0(path, "&ids=", ids)
  }
  if(!is.null(production_office)){
    path <- paste0(path, "&production-office=", production_office)
  }
  
  # Pagination
  if(!is.null(page)){
    path <- paste0(path, "&page=", page)
  }

  # Fields
  if(!is.null(fields)){
    path <- paste0(path, "&show-fields=", merge_multis(fields))
  }

  retrieved_data <- guardian_query(path, ...)[[1]]
  if(collapse){
    retrieved_data$results <- collapse_content(retrieved_data)
  }
  return(retrieved_data)
}
