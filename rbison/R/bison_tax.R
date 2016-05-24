#' Search for and collect taxonomic name data from the USGS Bison API using solr.
#'
#' @export
#'
#' @param query Name to search for. Required.
#' @param method The field to query by. See description below for details.
#' @param exact Exact matching or not. See examples. Defaults to FALSE.
#' @param parsed If TRUE (default) creates data.frame of names data output. Otherwise,
#' a list.
#' @param callopts Further args passed on to httr::GET for HTTP debugging/inspecting. In
#' \code{bison}, \code{bison_providers}, and \code{bison_stats}, \code{...} is used instead of
#' callopts, but \code{...} is used here to pass additional Solr params.
#' @param ... Further solr arguments passed in to the query. See examples below.
#' @return A list.
#' @description
#' See the SOLR documentation here \url{http://lucene.apache.org/solr/} for other
#' parameters you can use.
#'
#' The following two methods are possible, as far as I know you can only use one at
#' a time:
#' \itemize{
#'   \item vernacularName	The species specific common names that is searchable in a case
#'   insensitive way.
#'   \item scientificName	The species scientific name that is associated with a common
#'   name that is searchable in a case insensitive way.
#' }
#' @seealso \code{\link{bison_solr}} \code{\link{bison}}
#' @examples \dontrun{
#' # All taxa
#' bison_tax("*:*")
#' 
#' # Some example calls
#' bison_tax(query="*bear")
#' bison_tax(query="Helianthus", method="scientificName")
#'
#' # Exact argument, here nothing found with latter call as '*bear' doesn't exist,
#' # which makes sense
#' bison_tax(query="*bear", exact=FALSE)
#' bison_tax(query="*bear", exact=TRUE)
#'
#' # Using solr arguments (not all Solr arguments work)
#' ## Return a certain number of rows
#' bison_tax(query="*bear", method="vernacularName", rows=3)
#' ## Return certain fields
#' bison_tax(query="*bear", method="vernacularName", fl='vernacularName')
#'
#' # Curl options
#' library("httr")
#' bison_tax(query='*dolphin', callopts=verbose())
#' }

bison_tax <- function(query, method='vernacularName', exact=FALSE, parsed=TRUE,
  callopts=list(), ...)
{
  method <- match.arg(method, choices = c('vernacularName','scientificName'))
  if (!length(method) == 1) {
  	stop("method can only be of length 1")
  }
  url <- sprintf('http://bison.usgs.ornl.gov/solrstaging/%s/select', method)
  if (exact) { 
    qu_ <- paste0('"', query, '"') 
  } else { 
    qu_ <- query 
  }
  args <- bs_compact(list(q = qu_, wt = "json", ...))
  tt <- GET(url, query = args, c(config(followlocation = 1), callopts))
  stop_for_status(tt)
  out <- content(tt)
  temp <- list(
    numFound = out$response$numFound,
    names = out$response$docs,
    highlight = out$highlighting,
    facets = out$facet_counts
  )

  if (parsed) {
    data <- rbind_all(lapply(out$response$docs, data.frame, stringsAsFactors = FALSE))
    data$X_version_ <- NULL
    temp$names <- data
  }

  return( temp )
}
