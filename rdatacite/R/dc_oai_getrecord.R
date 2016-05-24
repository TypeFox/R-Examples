#' Retrieve an individual record from the DataCite repository
#'
#' @export
#'
#' @inheritParams dc_oai_listmetadataformats
#' @examples \dontrun{
#' # you can pass in just the numeric part of the ID
#' dc_oai_getrecord(56225)
#' dc_oai_getrecord(c(56225, 6667400))
#'
#' # Or, the entire thing
#' dc_oai_getrecord("oai:oai.datacite.org:56225")
#'
#' # Or, mixed
#' dc_oai_getrecord(c("56225", "oai:oai.datacite.org:6667400"))
#'
#' today <- format(Sys.Date(), "%Y-%m-%d")
#' temp <- dc_oai_listidentifiers(from = today)
#' dc_oai_getrecord(temp$identifier[1:2])
#' }

dc_oai_getrecord <- function(id) {
  tt <- sapply(id, function(z) grepl('oai:oai.datacite.org', z))
  if (!all(tt)) {
    id <- unname(sapply(id, strextract, pattern = "[0-9]+"))
    id <- paste("oai:oai.datacite.org:", id, sep = "")
  }
  oai::get_records(ids = id, prefix = "oai_dc", url = dc_oai_base())
}
