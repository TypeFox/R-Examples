#' List available metadata formats from various providers.
#'
#' @export
#' @template url_ddd
#' @param id The OAI-PMH identifier for the record. Optional.
#' @examples \dontrun{
#' list_metadataformats()
#'
#' # no metadatformats for an identifier
#' list_metadataformats(id = "oai:oai.datacite.org:22")
#'
#' # metadatformats available for an identifier
#' list_metadataformats(id = "oai:oai.datacite.org:32348")
#'
#' # curl options
#' library("httr")
#' list_metadataformats(id = "oai:oai.datacite.org:32348", config = verbose())
#' }
list_metadataformats <- function(url = "http://oai.datacite.org/oai", id = NULL, ...) {
  check_url(url)
  if (!is.null(id)) {
    setNames(lapply(id, one_mf, url = url, ...), id)
  } else {
    one_mf(id, url, ...)
  }
}

one_mf <- function(identifier, url, ...) {
  args <- sc(list(verb = 'ListMetadataFormats', identifier = identifier))
  res <- GET(url, query = args, ...)
  stop_for_status(res)
  out <- content(res, "text", encoding = "UTF-8")
  xml <- xml2::read_xml(out)
  rbind.fill(lapply(xml_children(xml_children(xml)[[3]]), get_headers))
}
