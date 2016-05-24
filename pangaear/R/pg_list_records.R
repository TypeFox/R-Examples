#' List records from Pangaea
#'
#' @export
#' @param prefix A character string to specify the metadata format in OAI-PMH requests issued to
#' the repository. The default (\code{"oai_dc"}) corresponds to the mandatory OAI unqualified
#' Dublin Core metadata schema.
#' @param from Character string giving datestamp to be used as lower bound for datestamp-based
#' selective harvesting (i.e., only harvest records with datestamps in the given range). Dates
#' and times must be encoded using ISO 8601. The trailing Z must be used when including time.
#' OAI-PMH implies UTC for data/time specifications.
#' @param until Character string giving a datestamp to be used as an upper bound,
#' for datestamp-based selective harvesting (i.e., only harvest records with datestamps in
#' the given range).
#' @param set A character string giving a set to be used for selective harvesting (i.e., only
#' harvest records in the given set).
#' @param token	(character) a token previously provided by the server to resume a
#' request where it last left off. 50 is max number of records returned. We will
#' loop for you internally to get all the records you asked for.
#' @param as (character) What to return. One of "df" (for data.frame; default),
#' "list", or "raw" (raw text)
#' @param ... Curl debugging options passed on to \code{\link[httr]{GET}}
#' @return XML character string, data.frame, or list, depending on what requested
#' witht the \code{as} parameter
#' @examples \dontrun{
#' pg_list_records(from='2015-09-01', until='2015-09-10')
#' pg_list_records(set='geomound', from='2015-01-01', until='2015-01-05')
#'
#' # When no results found > "'noRecordsMatch'"
#' # pg_list_records(set='geomound', from='2015-01-01', until='2015-01-01')
#'
#' # More examples
#' pg_list_records(set='citable', from='2015-10-01', until='2015-10-02')
#'
#' pg_list_records(prefix="iso19139", set='citable', from='2015-10-01', until='2015-10-02')
#' ## FIXME - below are broken
#' # pg_list_records(prefix="dif", set='citable', from='2015-01-01', until='2015-01-05')
#' # pg_list_records(prefix="dif", set='project4094', from='2015-01-01', until='2015-01-05')
#' }

pg_list_records <- function(prefix = "oai_dc", from = NULL, until = NULL, set = NULL,
                            token = NULL, as = "df", ...) {
  oai::list_records(url = baseoai(), prefix = prefix, from = from,
                    until = until, set = set, token = token, as = as, ...)
}
