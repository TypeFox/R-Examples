#' Use the DataCite RESTful API for search, facet, stats, and mlt queries.
#'
#' @export
#' @template search
#'
#' @examples \dontrun{
#' # Search
#' ## various searches
#' ### containing laser
#' dc_search(q = "laser", fl=c('doi','publicationYear'), rows=5)
#'
#' ### containing laser and safety
#' dc_search(q = "laser safety", fl=c('doi','publicationYear'), rows=5)
#'
#' ### containing laser and not medicine
#' dc_search(q = "laser -medicine", fl=c('doi','publicationYear'), rows=5)
#'
#' ### containing laser and red or green
#' dc_search(q = "laser AND (red OR green)", fl=c('doi','publicationYear'), rows=5)
#'
#' ### containing wind turbine as a phrase in the title
#' dc_search(q = 'title:"wind turbine"', fl=c('doi','title'), rows=5)
#'
#' ### containing geology in metadata field subject
#' dc_search(q = "subject:geology", fl=c('doi','subject'), rows=5)
#'
#' ### belonging to DataCite members BL or TIB
#' dc_search(q = "allocator:(BL OR TIB)", fl=c('doi','allocator'), rows=5)
#'
#' ### published between 2000 and 2005
#' dc_search(q = "publicationYear:[2000 TO 2005]", fl=c('doi','publicationYear'), rows=5)
#'
#' ### uploaded to DataCite in the last 5 days
#' dc_search(q = "uploaded:[NOW-5DAYS TO NOW]", fl=c('doi','uploaded'), rows=5)
#'
#' ## search with csv output
#' ### If you want to use wt=csv, install development version of solr package, see Details
#' dc_search(q = 'wind', fl=c('doi','title'), wt='csv')
#'
#' # Faceting
#' dc_facet(q = "wind", facet.field='publisher_facet')
#'
#' # Stats
#' dc_stats(q = "ecology", stats.field='date')
#'
#' # More like this (aka mlt)
#' dc_mlt(q = "ecology", mlt.fl='title', mlt.count=5, fl=c('doi','title'))
#' }
dc_search <- function(..., verbose = TRUE, errors = "simple", proxy = NULL, callopts=list()) {
  check_conn(verbose, errors, proxy)
  solrium::solr_search(..., callopts = callopts)
}

#' @export
#' @rdname dc_search
dc_facet <- function(..., verbose = TRUE, errors = "simple", proxy = NULL, callopts=list()) {
  check_conn(verbose, errors, proxy)
  solrium::solr_facet(..., callopts = callopts)
}

#' @export
#' @rdname dc_search
dc_stats <- function(..., verbose = TRUE, errors = "simple", proxy = NULL, callopts=list()) {
  check_conn(verbose, errors, proxy)
  solrium::solr_stats(..., callopts = callopts)
}

#' @export
#' @rdname dc_search
dc_mlt <- function(..., verbose = TRUE, errors = "simple", proxy = NULL, callopts=list()) {
  check_conn(verbose, errors, proxy)
  solrium::solr_mlt(..., callopts = callopts)
}

check_conn <- function(verbose, errors, proxy) {
  solrium::solr_connect(dc_base(), proxy = proxy, errors = errors, verbose = verbose)
}
