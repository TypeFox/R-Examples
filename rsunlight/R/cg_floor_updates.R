#' Search for data on upcoming bills.
#'
#' @export
#'
#' @param chamber The chamber this update took place in. 'house' or 'senate'.
#' @param timestamp The time this update took place. For Senate updates, this actually means the
#' time our system first observed the update, and is susceptible to error; the Senate does not
#' offer precise timestamps.
#' @param congress The number of the Congress this update took place during.
#' @param legislative_day The 'legislative day' this update took place in. The 'legislative day'
#' is a formal construct that is usually, but not always, the same as the calendar day. For
#' example, if a day's session of Congress runs past midnight, the legislative_day will often
#' stay the same as it was before midnight, until that session adjourns. On January 3rd, it is
#' possible that the same legislative_day could span two Congresses. (This occurred in 2013.)
#' @param year The 'legislative year' of the update. This is not quite the same as the calendar
#' year - the legislative year changes at noon EST on January 3rd. A vote taken on January 1, 2013
#' has a 'legislative year' of 2012.
#' @param bill_ids An array of IDs of bills that are referenced by or associated with this floor
#' update.
#' @param roll_ids An array of IDs of roll call votes that are referenced by or associated with
#' this floor update.
#' @param legislator_ids An array of bioguide IDs of legislators that are referenced by this floor
#' update.
#'
#' @template cg
#' @template cg_query
#' @examples \dontrun{
#' cg_floor_updates()
#' cg_floor_updates(chamber='house', query='Agreed to by voice vote')
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_floor_updates(chamber = c("house", "senate"))
#' }

cg_floor_updates <- function(chamber=NULL, timestamp=NULL, congress=NULL, legislative_day=NULL,
  year=NULL, bill_ids=NULL, roll_ids=NULL, legislator_ids=NULL, query=NULL,
  fields=NULL, page=1, per_page=20, order=NULL, key = NULL, as = 'table', ...) {

  key <- check_key(key)
  args <- sc(list(apikey=key, chamber=chamber, timestamp=timestamp, congress=congress,
                          legislative_day=legislative_day, year=year, bill_ids=bill_ids,
                          roll_ids=roll_ids, legislator_ids=legislator_ids, query=query,
                          per_page=per_page, page=page, fields=fields, order=order))
  give_cg(as, cgurl(), "/floor_updates", args, ...)
}
