#' Get basic title, item, and page metadata for each page on which the specified
#' name appears.
#'
#' @export
#' @param namebankid (not used if 'name' specified) NameBank identifier for a
#' name (numeric)
#' @param name (not used if 'namebankid' specified) a complete name string (character)
#' @inheritParams bhl_getcollections
#'
#' @examples \dontrun{
#' bhl_namegetdetail(namebankid = 3501464)
#' bhl_namegetdetail(name = 'poa annua')
#' }

bhl_namegetdetail <- function(namebankid = NULL, name = NULL, as='table', key = NULL, ...) {
  args <- bhlc(list(op = "NameGetDetail", apikey = check_key(key),
                       namebankid = namebankid, name = name, format = as_f(as)))
  bhl_GET(as, args, ...)
}
