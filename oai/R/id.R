#' Identify the OAI-PMH service for each data provider.
#'
#' @export
#' @template url_ddd
#' @examples \dontrun{
#' # datacite
#' id("http://oai.datacite.org/oai")
#'
#' # arxiv
#' id("http://export.arxiv.org/oai2")
#'
#' # GBIF - http://www.gbif.org/
#' id("http://api.gbif.org/v1/oai-pmh/registry")
#'
#' # curl options
#' library("httr")
#' id("http://export.arxiv.org/oai2", config = verbose())
#' }
id <- function(url, ...) {
  check_url(url)
  rbind.fill(lapply(url, id_, ...))
}

id_ <- function(x, ...) {
  res <- GET(x, query = list(verb = "Identify"), ...)
  tt <- content(res, as = "text", encoding = "UTF-8")
  get_headers(xml_children(xml2::read_xml(tt))[[3]])
}
