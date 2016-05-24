#' Get records
#'
#' @export
#'
#' @template url_ddd
#' @template as
#' @param ids The OAI-PMH identifier for the record. One or more. Required.
#' @param prefix specifies the metadata format that the records will be
#'     returned in. Default: oai_dc
#' @examples \dontrun{
#' get_records("oai:oai.datacite.org:32255")
#' get_records(c("oai:oai.datacite.org:32255", "oai:oai.datacite.org:32325"))
#'
#' # Get a list
#' get_records("oai:oai.datacite.org:32255", as = "list")
#'
#' # Get raw text
#' get_records("oai:oai.datacite.org:32255", as = "raw")
#'
#' # from arxiv.org
#' get_records("oai:arXiv.org:0704.0001", url = "http://export.arxiv.org/oai2")
#'
#' # GBIF - http://www.gbif.org/
#' get_records(c("816f4734-6b49-41ab-8a1d-1b21e6b5486d", "95e3042f-f48d-4a04-8251-f755bebeced6"),
#'    url = "http://api.gbif.org/v1/oai-pmh/registry")
#' }
get_records <- function(ids, prefix = "oai_dc", url = "http://oai.datacite.org/oai", as = "df", ...) {
  check_url(url)
  out <- lapply(ids, each_record, url = url, prefix = prefix, as = as, ...)
  oai_give(do.call("c", out), as, "GetRecord")
}

each_record <- function(identifier, url, prefix, as, ...) {
  args <- sc(list(verb = "GetRecord", metadataPrefix = prefix, identifier = identifier))
  res <- GET(url, query = args, ...)
  stop_for_status(res)
  tt <- content(res, "text", encoding = "UTF-8")
  xml_orig <- xml2::read_xml(tt)
  handle_errors(xml_orig)
  if (as == "raw") {
    tt
  } else {
    xml <- xml2::xml_children(xml2::xml_children(xml_orig)[[3]])
    get_data(xml, as = as)
  }
}
