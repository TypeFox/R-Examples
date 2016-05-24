#' List OAI-PMH identifiers
#'
#' @export
#' @template url_ddd
#' @template as
#' @param prefix Specifies the metadata format that the records will be
#'     returned in.
#' @param from specifies that records returned must have been created/update/deleted
#'     on or after this date.
#' @param until specifies that records returned must have been created/update/deleted
#'     on or before this date.
#' @param set specifies the set that returned records must belong to.
#' @param token a token previously provided by the server to resume a request
#'     where it last left off.
#' @examples \dontrun{
#' # from
#' today <- format(Sys.Date(), "%Y-%m-%d")
#' list_identifiers(from = today)
#'
#' # from and until
#' list_identifiers(from = '2011-06-01T', until = '2011-07-01T')
#'
#' # longer time span
#' list_identifiers(from = '2011-06-01T', until = '2011-09-01T')
#'
#' # set parameter - here, using ANDS - Australian National Data Service
#' list_identifiers(from = '2011-09-01T', until = '2012-09-01T', set = "ANDS")
#'
#' # Get a list
#' list_identifiers(from = today, as = "list")
#'
#' # Get raw text
#' list_identifiers(from = today, as = "raw")
#' }
list_identifiers <- function(url = "http://oai.datacite.org/oai", prefix = "oai_dc", from = NULL,
                             until = NULL, set = NULL, token = NULL, as = "df", ...) {
  check_url(url)
  if (!is.null(token)) from <- until <- set <- prefix <- NULL
  args <- sc(list(verb = "ListIdentifiers", metadataPrefix = prefix, from = from,
                  until = until, set = set, resumptionToken = token))
  out <- while_oai(url, args, token, as, ...)
  oai_give(out, as, "ListRecords")
}
