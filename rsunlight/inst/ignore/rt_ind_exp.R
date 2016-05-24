#' Real-time Federal Campaign Finance API - New Filings
#'
#' @export
#' @param support_oppose_checked Whether the expenditures supports (S) or opposes (O)
#' a candidate.
#' @param min_spent The minimum money spent on an independent expenditure.
#' @param filer_committee_id_number The FEC id of the committee making the expenditure.
#' This is a nine-character code that begins with 'C'.
#' @param candidate_id_checked The FEC id of the candidate targeted (either supported or
#' opposed) by the independent expenditure. This is a nine-character code that begins
#' with 'H', 'S', or 'P'.
#' @param district_checked The internal id of the district the candidate is running in.
#' The id can be determined by using the /districts/ endpoint.
#' @param candidate_state_checked The two-digit postal code of the state in which the
#' candidate targeted is running, e.g. 'NC' or 'OK'
#' @param candidate_office_checked A one-digit code for the office the candidate targeted
#' by the expenditure is seeking: 'H' for House or 'S' for Senate.
#' @param candidate_district_checked A two digit code showing the district number of
#' house candidates only. The first district would be '01'. States with at-large house
#' districts, e.g. North Dakota, are assigned '00'. This number is not reliably assigned
#' to senate districts; sometimes a '00' is used, sometimes the value is null
#' @param candidate_party_checked A one digit code representing the party of the candidate
#' targeted by the expenditure: 'D' for Democrat and 'R' for Republican. This value is
#' often missing for third party candidates.
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
#' rt_ind_exp(candidate_party_checked = 'D', page_size = 2)
#' rt_ind_exp(min_spent = 10000000)
#' }
rt_ind_exp <- function(support_oppose_checked = NULL, min_spent = NULL,
  filer_committee_id_number = NULL, candidate_id_checked = NULL, district_checked = NULL,
  candidate_state_checked = NULL, candidate_office_checked = NULL,
  candidate_district_checked = NULL, candidate_party_checked = NULL, as = 'table',
  page = 1, page_size = 10, key = NULL, ...) {

  key <- check_key(key)
  args <- sc(list(apikey = key, format = "json", support_oppose_checked = support_oppose_checked,
                  min_spent = min_spent, filer_committee_id_number = filer_committee_id_number,
                  candidate_id_checked = candidate_id_checked, district_checked = district_checked,
                  candidate_state_checked = candidate_state_checked,
                  candidate_office_checked = candidate_office_checked,
                  candidate_district_checked = candidate_district_checked,
                  candidate_party_checked = candidate_party_checked,
                  page = page, page_size = page_size))
  return_obj(as, query(paste0(rtieurl(), "/independent-expenditures/"), args, ...))
}
