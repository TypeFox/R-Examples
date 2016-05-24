#' Search the Dryad Solr endpoint.
#'
#' @export
#'
#' @param ... Parameters to pass to any of the \code{solr_*} functions in the \code{\link{solr}}
#' package.
#' @examples \dontrun{
#' # Basic search
#' d_solr_search(q="Galliard")
#'
#' # Basic search, restricting to certain fields
#' d_solr_search(q="Galliard", fl='handle,dc.title_sort')
#'
#' # Search all text for a string, but limits results to two specified fields:
#' d_solr_search(q="dwc.ScientificName:drosophila", fl='handle,dc.title_sort')
#'
#' # Dryad data based on an article DOI:
#' d_solr_search(q="dc.relation.isreferencedby:10.1038/nature04863",
#'    fl="dc.identifier,dc.title_ac")
#'
#' # All terms in the dc.subject facet, along with their frequencies:
#' d_solr_facet(q="location:l2", facet.field="dc.subject_filter", facet.minCount=1,
#'    facet.limit=10)
#'
#' # Article DOIs associated with all data published in Dryad over the past 90 days:
#' d_solr_search(q="dc.date.available_dt:[NOW-90DAY/DAY TO NOW]",
#'    fl="dc.relation.isreferencedby", rows=10)
#'
#' # Data DOIs published in Dryad during January 2011, with results returned in JSON format:
#' query <- "location:l2 dc.date.available_dt:[2011-01-01T00:00:00Z TO 2011-01-31T23:59:59Z]"
#' d_solr_search(q=query, fl="dc.identifier", rows=10)
#'
#' # Highlight
#' d_solr_highlight(q="bird", hl.fl="dc.description")
#'
#' # More like this
#' d_solr_mlt(q="bird", mlt.count=10, mlt.fl='dc.title_sort', fl='handle,dc.title_sort')
#'
#' # Stats
#' d_solr_stats(q="*:*", stats.field="dc.date.accessioned.year")
#' }
d_solr_search <- function(...){
  solr_search(..., base = dsolrbase())
}

#' @export
#' @rdname d_solr_search
d_solr_facet <- function(...){
  solr_facet(..., base = dsolrbase())
}

#' @export
#' @rdname d_solr_search
d_solr_group <- function(...){
  solr_group(..., base = dsolrbase())
}

#' @export
#' @rdname d_solr_search
d_solr_highlight <- function(...){
  solr_highlight(..., base = dsolrbase())
}

#' @export
#' @rdname d_solr_search
d_solr_mlt <- function(...){
  solr_mlt(..., base = dsolrbase())
}

#' @export
#' @rdname d_solr_search
d_solr_stats <- function(...){
  solr_stats(..., base = dsolrbase())
}

dsolrbase <- function() "http://datadryad.org/solr/search/select"
