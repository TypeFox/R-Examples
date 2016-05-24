#' Search Legislators on OpenStates.
#'
#' @param state state two-letter abbreviation (character)
#' @param first_name first name of legislator (character)
#' @param last_name last name of legislator (character)
#' @param chamber one of 'upper' or 'lower' (character)
#' @param active TRUE or FALSE (character)
#' @param term filter by legislators who served during a certain term (character)
#' @param district legislative district (character)
#' @param party democratic or republican (character)
#' @param fields You can request specific fields by supplying a vector of fields names. Many fields
#' are not returned unless requested. If you don't supply a fields parameter, you will get the
#' most commonly used subset of fields only. To save on bandwidth, parsing time, and confusion,
#' it's recommended to always specify which fields you will be using.
#' @param as (character) One of table (default), list, or response (httr response object)
#' @param key your SunlightLabs API key; loads from .Rprofile
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @return List of output fields.
#' @export
#' @examples \dontrun{
#' os_legislatorsearch(state = 'ca', party = 'democratic')
#' os_legislatorsearch(state = 'tx', party = 'democratic', active = TRUE)
#' os_legislatorsearch(state = 'nv', party = 'republican')
#' os_legislatorsearch(state = 'dc', chamber = 'upper')
#'
#' # pass in more than one value for some parameters
#' os_legislatorsearch(state = c('dc', 'or'), chamber = 'upper')
#' os_legislatorsearch(first_name = c('jane', 'bob'), chamber = 'upper')
#' }
os_legislatorsearch <- function(state = NULL, first_name = NULL,
    last_name = NULL, chamber = NULL, active = NULL, term = NULL, district = NULL,
    party = NULL, fields = NULL, as = 'table', key = NULL, ...) {

  key <- check_key(key)
  args <- sc(list(apikey = key, state = state, first_name = first_name,
                       last_name = last_name, chamber = chamber, active = active,
                       term = term, district = district, party = party,
                       fields = paste(fields, collapse = ",")))
  give(as, osurl(), "/legislators", args, ...)
}
