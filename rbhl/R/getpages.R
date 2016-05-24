#' Get many OCR-generated pages given a single item id
#'
#' @export
#' @param itemid the item id (character). Required
#' @param key Your BHL API key, either enter, or loads from your \code{.Renviron} as \code{BHL_KEY}
#' or from \code{.Rprofile} as \code{bhl_key}.
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @examples \dontrun{
#' books <- bhl_booksearch(title='Selborne', lname='White', volume=2, edition='new', year=1825,
#'    collectionid=4, language='eng')
#' getpages(itemid=16800)
#' }

getpages <- function(itemid, key = NULL, ...){
  res <- bhl_getitempages(itemid, ...)
  out <- lapply(res$data$PageID, function(x){
    tmp <- bhl_getpageocrtext(page = x)
    if (grepl("OCR text unavailable for this page", tmp))
      NULL
    else
      tmp
  })
  setNames(out, res$PageID)
}
