#' Return metadata about a page.
#'
#' You may choose to include the OCR text and a list of names found on the page.
#'
#' @export
#' @param page page number to get
#' @param ocr return ocr text of the page (TRUE/FALSE)
#' @param names return the names that appear on the page (TRUE/FALSE)
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getpagemetadata(page=1328690, ocr=TRUE)
#' bhl_getpagemetadata(page=1328690, ocr=TRUE, as='json')
#' bhl_getpagemetadata(page=1328690, ocr=TRUE, as='xml')
#' bhl_getpagemetadata(page=1328690, ocr=TRUE, as='list')
#' }

bhl_getpagemetadata <- function(page = NULL, ocr = FALSE, names = FALSE,
  as = 'table', key = NULL, ...)
{
  args <- bhlc(list(op = "GetPageMetadata", apikey = check_key(key), pageid=page,
                       ocr=if(ocr) 't' else NULL, names=if(names) 't' else NULL,
                       format=as_f(as)))
  bhl_GET(as, args, ...)
}
