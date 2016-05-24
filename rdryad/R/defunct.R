#' Download metadata for all Dryad oai's for defined time period
#'
#' This function is defunct
#'
#' @keywords internal
getalldryad_metadata <- function(...) {
  .Defunct(msg = "This function is defunct. Use OAI-PMH via dr_*() functions or Solr based search via d_*() functions")
}

#' Search metadata for search terms using regex
#'
#' This function is defunct
#'
#' @keywords internal
search_dryad <- function(...) {
  .Defunct(msg = "This function is defunct. Use OAI-PMH via dr_*() functions or Solr based search via d_*() functions")
}

#' Download metadata for individual Dryad id's
#'
#' This function changed name to \code{\link{dr_get_records}}
#'
#' @keywords internal
download_dryadmetadata <- function(...) {
  .Defunct(msg = "This function is defunct. Use dr_get_records()", new = "dr_get_records", package = "rdryad")
}
