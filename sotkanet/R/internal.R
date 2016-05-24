
#' Description:
#' SotkanetCollect converts the list object from Sotkanet to a data.frame
#'
#' Arguments:
#'   @param x input data (from SotkanetIndicators or SotkanetRegions etc.)
#'   @param name name for the column ("indicator", "region", etc.)
#'
#' Returns:
#'   @return sotkanet data table
#'
#' @references
#' See citation("sotkanet") 
#' @author Einari Happonen / Opasnet. Maintainer: Louhos \email{louhos@@googlegroups.com}
#' @keywords utilities

SotkanetCollect <- function(x, name) {

  if (length(x$id) == 1) {
    x <- list(x)
  }

  if (name == "region") {
    out <- data.frame(list(
      region = sapply(x, function (xi) {xi$id}),
      region.title.fi = sapply(x, function (xi) {xi$title[["fi"]]}),
      region.code = sapply(x, function (xi) {xi$code}),
      region.category = sapply(x, function (xi) {xi$category}),
      region.uri = gsub("NULL", "", sapply(x, function (xi) {xi$uri}))))
  } else if (name == "indicator") {
    out <- data.frame(list(
      indicator = sapply(x, function (xi) {xi$id}),
      indicator.title.fi = sapply(x, function (xi) {xi$title[["fi"]]}),
      indicator.organization = sapply(x, function (xi) {xi$organization$id}),
      indicator.organization.title.fi = sapply(x, function (xi) {xi$organization$title[["fi"]]})
      #indicator.last.update = gsub("NULL", "", sapply(x, function (xi) {xi[["data-update"]]}))
    ))
  }

  return(out)
}

