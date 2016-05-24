#' Search OpenStates bills.
#'
#' @param terms search terms bill search (character)
#' @param state state two-letter abbreviation (character)
#' @param window a string representing what time period to search across.
#'    Pass 'session' to search bills from the state's current or most
#'    recent legislative session, 'term' to search the current or most
#'    recent term, 'all' to search as far back as the Open State Project
#'    has data for, or supply 'session:SESSION_NAME' or 'term:TERM_NAME'
#'    (e.g. 'session:2009' or 'term:2009-2010') to search a specific
#'    session or term.
#' @param chamber one of 'upper' or 'lower' (character)
#' @param sponsor_id only return bills sponsored by the legislator with the
#'    given id (corresponds to leg_id)
#' @param updated_since only return bills that have been updated since a given date,
#'    YYYY-MM-DD format
#' @param subject filter by bills that are about a given subject. If multiple
#'    subject parameters are supplied then only bills that match all of
#'    them will be returned. See list of subjects
#' @param type (character) Only bills of a given type (e.g. 'bill', 'resolution', etc.)
#' @param search_window By default all bills are searched, but if a time window is desired the
#' following options can be passed to search_window: all (Default, include all sessions.),
#' term (Only bills from sessions within the current term.), session (Only bills from the current
#' session.), session:2009 (Only bills from the session named 2009.), term:2009-2011 (Only bills
#' from the sessions in the 2009-2011 session.)
#' @param sort (character) One of 'first' (default), 'last', 'signed', 'passed_lower',
#'    'passed_upper', 'updated_at', or 'created_at'.
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
#' @return a data.frame of bills.
#' @export
#' @examples \dontrun{
#' os_billsearch(terms = 'agriculture', state = 'tx')
#' os_billsearch(terms = 'agriculture', state = 'tx', chamber = 'upper')
#' os_billsearch(terms = 'taxi', state = 'dc')
#' os_billsearch(terms = 'taxi', state = 'dc', per_page=3)
#' os_billsearch(terms = 'taxi', state = 'dc', per_page=3, sort='created_at')
#' os_billsearch(terms = 'taxi', state = 'dc', type='resolution')
#'
#' # Search window
#' os_billsearch(terms = 'climate change', search_window='term')
#' os_billsearch(terms = 'climate change', search_window='term:2009-2011')
#' os_billsearch(terms = 'climate change', search_window='session')
#' os_billsearch(terms = 'climate change', search_window='session:2009')
#'
#' os_billsearch(terms = 'agriculture', state = 'tx', per_page=2)
#' os_billsearch(terms = 'agriculture', state = 'tx', per_page=2, page=2)
#' os_billsearch(terms = 'agriculture', state = 'tx', fields=c('id','created_at'), per_page=10)
#'
#' # Pass in more than one value for some parameters
#' os_billsearch(terms = 'agriculture', state = c('tx', 'or'))
#' os_billsearch(terms = 'agriculture', state = "or", chamber = c('upper', 'lower'))
#' }
os_billsearch <- function(terms = NULL, state = NULL, window = NULL,
    chamber = 'upper', sponsor_id = NULL, updated_since = NULL, subject = NULL, type=NULL,
    search_window=NULL, sort=NULL, page=NULL, per_page=NULL, fields = NULL, as ='table',
    key = NULL, ...) {

  key <- check_key(key)
  args <- sc(list(apikey = key, q = terms, state = state, window = window,
                       chamber = chamber, sponsor_id = sponsor_id,
                       updated_since = updated_since, subject = subject, type = type,
                       search_window = search_window, sort = sort, page = page, per_page = per_page,
                       fields = paste(fields, collapse = ",")))
  give(as, osurl(), "/bills", args, ...)
}
