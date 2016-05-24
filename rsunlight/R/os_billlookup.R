#' Lookup bills on OpenStates.
#'
#' @importFrom stringr str_extract ignore.case
#'
#' @param state state two-letter abbreviation (character), required
#' @param session session of congress (integer), e.g., 2009-2010 = 20092010,
#'    required
#' @param bill_id One or more identification numbers of bills (character), required
#' @param fields You can request specific fields by supplying a vector of fields names. Many fields
#' are not returned unless requested. If you don't supply a fields parameter, you will get the
#' most commonly used subset of fields only. To save on bandwidth, parsing time, and confusion,
#' it's recommended to always specify which fields you will be using.
#' @param per_page Number of records to return. Default: 20. Max: 50.
#' @param page Page to return. Default: 1. You can use this in combination with the
#' per_page parameter to get more than the default or max number of results per page.
#' @param as (character) One of table (default), list, or response (httr response object)
#' @param key your SunlightLabs API key; loads from .Rprofile
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @return List of output fields.
#' @export
#' @examples \dontrun{
#' os_billlookup(state='ca', session=20092010, bill_id='AB 667')
#' os_billlookup(state='ca', session=20092010, bill_id='AB 667', per_page=1)
#' os_billlookup(state='ca', session=20092010, bill_id='AB 667', per_page=1, fields='id')
#' os_billlookup(state='ca', session=20092010, bill_id='AB 667',
#'    per_page=3, fields=c('id','title'))
#' os_billlookup(state='ca', session=20092010, bill_id='SB 425')
#' os_billlookup(state='ca', session=20092010, bill_id=c('AB 667','SB 425'))
#'
#' library('httr')
#' os_billlookup(state='ca', session=20092010, bill_id='AB 667', config=verbose(), per_page=1)
#' }

os_billlookup <- function(state = NULL, session = NULL, bill_id = NULL,
  fields = NULL, per_page = NULL, page = NULL, as = 'table', key = NULL, ...) {

  key <- check_key(key)
  bills <- NULL
  if (length(bill_id) > 1) {
    bills <- paste(bill_id, collapse = "|")
    bill_id <- NULL
  }
  args <- sc(list(apikey = key, state = state, session = session, bill_id = bill_id,
          bill_id__in = bills, per_page = per_page, page = page, fields = paste(fields, collapse = ",")))
  return_obj(as, query(paste0(osurl(), "/bills"), args, ...))
}
