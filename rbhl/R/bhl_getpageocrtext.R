#' Return the OCR-generated text of a page.
#'
#' @export
#' @param page page number to get
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getpageocrtext(page=1328690)
#' bhl_getpageocrtext(1328690, 'json')
#' bhl_getpageocrtext(1328690, 'xml')
#' }

bhl_getpageocrtext <- function(page = NULL, as = 'table', key = NULL, ...)
{
  args <- bhlc(list(op = "GetPageOcrText", apikey = check_key(key), format=as_f(as), pageid=page))
  bhl_GET(as, args, ...)
}
