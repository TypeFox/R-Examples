#' Search for data on upcoming bills.
#'
#' @export
#'
#' @param scheduled_at The exact time at which our systems first spotted this bill
#' on the schedule in this chamber and on this legislative day. Currently, we check
#' the schedules every 15 minutes.
#' @param legislative_day The legislative day this bill is scheduled for. Combine
#' with the range field to understand precision. May be null.
#' @param range How precise this information is. One of 'day', 'week', or NULL. See
#' Details for more.
#' @param congress The number of the Congress this bill has been scheduled in.
#' @param chamber The chamber which has scheduled this bill.
#' @param source_type The source for this information. 'house_daily' (Majority Leader daily
#' schedule or 'senate_daily' (Senate Democrats' Floor feed.
#' @param bill_id The ID of the bill that is being scheduled.
#'
#' @template cg
#' @template cg_query
#'
#' @details More info for range parameter
#' \itemize{
#'  \item day: bill has been scheduled specifically for the legislative_day.
#'  \item week: bill has been scheduled for the 'Week of' the legislative_day.
#'  \item NULL: bill has been scheduled at an indefinite time in the future. (legislative_day
#'  is null.)
#' }
#'
#' The 'legislative day' is a formal construct that is usually, but not always, the same as the
#' calendar day. For example, if a day's session of Congress runs past midnight, the
#' legislative_day will often stay the same as it was before midnight, until that session adjourns.
#' On January 3rd, it is possible that the same legislative_day could span two Congresses.
#' (This occurred in 2013.)
#'
#' @examples \dontrun{
#' cg_upcoming_bills()
#'
#' # Pass in more than one value to a parameter
#' cg_upcoming_bills(chamber = c("house", "senate"))
#' }

cg_upcoming_bills <- function(scheduled_at=NULL, legislative_day=NULL, range=NULL, congress=NULL,
  chamber=NULL, source_type=NULL, bill_id=NULL, query=NULL, fields=NULL, page=1, per_page=20, order=NULL,
  key = NULL, as = 'table', ...) {

  key <- check_key(key)
  args <- sc(list(apikey=key, scheduled_at=scheduled_at, legislative_day=legislative_day,
                          range=range, congress=congress, chamber=chamber, source_type=source_type,
                          bill_id=bill_id, per_page=per_page, page=page, fields=fields,
                          query=query, order=order))
  give_cg(as, cgurl(), "/upcoming_bills", args, ...)
}
