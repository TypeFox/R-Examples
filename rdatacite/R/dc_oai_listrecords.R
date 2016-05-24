#' List the records in the DataCite metadata repository.
#'
#' @export
#' @inheritParams dc_oai_listidentifiers
#' @examples \dontrun{
#' dc_oai_listrecords(from = '2011-06-01T', until = '2011-07-01T')
#' dc_oai_listrecords(from = '2011-06-01T', until = '2011-11-01T', set = "REFQUALITY")
#' }
dc_oai_listrecords <- function(from = NULL, until = NULL, set = NULL, prefix = 'oai_dc',
                               token = NULL, ...) {
  oai::list_records(dc_oai_base(), prefix = prefix, from = from, until = until, set = set, ...)
}
