#' Real-time Federal Campaign Finance API - New Filings
#'
#' @export
#' @param candidate_id (character) The FEC's id for the candidate connected to this committee.
#' Only authorized committees are connected to candidates; leadership pacs are unconnected.
#' @param committee_class (character) The one letter code that describes the type of committee (see
#' the committee reference above). Multiple committee types are allowed--a search for
#' house and senate candidate committees would use HS, for instance.
#' @param min_coh (numeric) The minimum cash on hand of a report to include. Not all electronic
#' reports include a cash on hand, so using this option will exclude all of those that
#' do not report this figure.
#' @param min_raised (numeric) The minimum money raised. Not all electronic reports include an
#' amount raised, so using this option will exclude all of those that do not report this
#' figure.
#' @param min_spent (numeric) The minimum money spent. Not all electronic reports include an amount
#' spent, so using this option will exclude all of those that do not report this figure.
#' @param period_type (character) The type of periodic report to include--could be monthly, quarterly
#' or semiannual. The value includes the type; M for monthly, Q for quarterly, and S for
#' semiannually and the period number; so a second quarter report would be Q2 and a monthly
#' report covering August would be M8. Using this option will only return periodic reports.
#' @param year_covered (numeric) The calendar year that a periodic report covers. Non-periodic reports
#' that do not cover a particular time period will be excluded.
#' @param time_range (character) Only show reports received within one of the following time ranges:
#' today, in the last seven days, or at any time in this cycle.
#' @param report_type The type of report to include. See Details.
#' @param is_superceded (logical) When an amended version of an older filing is received
#' we consider the original to have been superseded. When set to false, show only the
#' most recent versions of an electronic filing; when set to True, show only original
#' filings which have since been amended.
#' @param format (character) The format of the returned data. json (default) or csv
#' @param page (integer) The page number of results to return. Default is 1.
#' @param page_size (integer) The number of page results to return. for json, the default is
#' 100 and the maximum is 100. This parameter is ignored for csv downloads, which
#' return all matching rows, with a limit of 2000 total rows returned.
#' @return Data frame of observations by date.
#' @details report_type options:
#' \itemize{
#'  \item monthly: Monthly/quarterly reports
#'  \item ies: Independent/coordinated expenditures
#'  \item F6:48-hour notice of contributions/loans received
#'  \item F9:24-hour notice of electioneering communications
#'  \item F2:Statement of candidacy
#'  \item F4:Convention committee report
#'  \item F13:Inaugural committee report
#' }
#' @examples \dontrun{
#' rt_new_filing(year_covered = 2014, page_size = 2)
#' rt_new_filing(year_covered = 2014, min_spent = 10000)
#' rt_new_filing(report_type = 'monthly')
#' }
rt_new_filing <- function(candidate_id = NULL, committee_class = NULL, min_coh = NULL,
  min_raised = NULL, min_spent = NULL, period_type = NULL, year_covered = NULL,
  time_range = NULL, report_type = NULL, is_superceded = NULL, as = 'table',
  page = 1, page_size = 10, key = NULL, ...) {

  key <- check_key(key)
  args <- sc(list(apikey = key, format = "json", candidate_id = candidate_id,
                  committee_class = committee_class, min_coh = min_coh,
                  min_raised = min_raised, min_spent = min_spent,
                  period_type = period_type, year_covered = year_covered,
                  time_range = time_range, report_type = report_type,
                  is_superceded = is_superceded, page = page,
                  page_size = page_size))
  return_obj(as, query(paste0(rtieurl(), "/new_filing/"), args, ...))
}
