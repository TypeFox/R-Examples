#' Gets details for bills.
#'
#' Data on bills in Congress goes back to 2009, and comes from a mix of sources:
#' \itemize{
#'  \item Scrapers at github.com/unitedstates for most data, including core status and
#'  history information.
#'  \item Bulk data at GPO's FDSys for version information, and full text.
#'  \item The House' MajorityLeader.gov and Senate Democrats' official site for
#'  notices of upcoming debate.
#' }
#'
#' @export
#' @template bills
#' @template cg
#'
#' @details
#' History: The history field includes useful flags and dates/times in a bill's life.
#' The above is a real-life example of H.R. 3590 - not all fields will be present
#' for every bill. Time fields can hold either dates or times - Congress is inconsistent
#' about providing specific timestamps.
#'
#' @return Committee details including subcommittees and all members.
#'
#' @examples \dontrun{
#' # Bill lookup (i.e., filter)
#' cg_bills(congress=113, history.enacted=TRUE)
#' cg_bills(history.active=TRUE, order='last_action_at')
#' cg_bills(sponsor.party='R', history.vetoed=TRUE)
#' cg_bills(enacted_as.law_type='private', order='history.enacted_at')
#' cg_bills(bill_type__in='hjres|sjres', history.house_passage_result__exists=TRUE,
#'    history.senate_passage_result__exists=TRUE)
#'
#' # Bill search
#' cg_bills(query='health care')
#' cg_bills(query='health care', history.enacted=TRUE)
#' cg_bills(query='freedom of information')
#' cg_bills(query='"freedom of information" accountab*')
#' cg_bills(query="'transparency accountability'~5", highlight=TRUE)
#'
#' # Disable pagination
#' cg_bills(per_page='all')
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_bills(bill_id = c("hjres131-113", "hjres129-113", "s2921-113"))
#' }

cg_bills <- function(query = NULL, bill_id = NULL, bill_type = NULL, number = NULL,
  congress = NULL, chamber = NULL, introduced_on = NULL, last_action_at = NULL,
  last_vote_at = NULL, last_version_on = NULL, highlight = NULL, history.active = NULL,
  history.active_at = NULL, history.house_passage_result = NULL,
  history.house_passage_result_at = NULL, history.senate_cloture_result = NULL,
  history.senate_cloture_result_at = NULL, history.senate_passage_result = NULL,
  history.senate_passage_result_at = NULL, history.vetoed = NULL, history.vetoed_at = NULL,
  history.house_override_result = NULL, history.house_override_result_at = NULL,
  history.senate_override_result = NULL, history.senate_override_result_at = NULL,
  history.awaiting_signature = NULL, history.awaiting_signature_since = NULL,
  history.enacted = NULL, history.enacted_at = NULL,
  sponsor.party = NULL, enacted_as.law_type = NULL, bill_type__in = NULL,
  history.house_passage_result__exists = NULL, history.senate_passage_result__exists = NULL,
  nicknames=NULL, keywords=NULL, sponsor_id=NULL, cosponsor_ids=NULL, cosponsors_count=NULL,
  withdrawn_cosponsors_count=NULL, withdrawn_cosponsor_ids=NULL, committee_ids=NULL,
  related_bill_ids=NULL, enacted_as.congress=NULL,
  enacted_as.number=NULL, fields=NULL, page = 1, per_page = 20, order = NULL,
  key = NULL, as = 'table', ...) {

  key <- check_key(key)
  if (is.null(query)) {
    url <- paste0(cgurl(), '/bills')
  } else {
    url <- paste0(cgurl(), '/bills/search')
  }
  args <- sc(list(apikey=key,query=query,bill_id=bill_id,bill_type=bill_type,
    number=number,congress=congress,chamber=chamber,introduced_on=introduced_on,
    last_action_at=last_action_at,last_vote_at=last_vote_at,last_version_on=last_version_on,
    highlight=highlight,history.active=ll(history.active), history.active_at=history.active_at,
    history.house_passage_result=history.house_passage_result,
    history.house_passage_result_at=history.house_passage_result_at,
    history.senate_cloture_result=history.senate_cloture_result,
    history.senate_cloture_result_at=history.senate_cloture_result_at,
    history.senate_passage_result=history.senate_passage_result,
    history.senate_passage_result_at=history.senate_passage_result_at,
    history.vetoed=ll(history.vetoed), history.vetoed_at=history.vetoed_at,
    history.house_override_result=history.house_override_result,
    history.house_override_result_at=history.house_override_result_at,
    history.senate_override_result=history.senate_override_result,
    history.senate_override_result_at=history.senate_override_result_at,
    history.awaiting_signature=ll(history.awaiting_signature),
    history.awaiting_signature_since=history.awaiting_signature_since,
    history.enacted=ll(history.enacted), history.enacted_at=history.enacted_at,
    sponsor.party=sponsor.party,order=order,enacted_as.law_type=enacted_as.law_type,
    bill_type__in=bill_type__in,
    history.house_passage_result__exists=history.house_passage_result__exists,
    history.senate_passage_result__exists=history.senate_passage_result__exists,
    page=page,per_page=per_page,fields=fields,
    nicknames=nicknames, keywords=keywords, sponsor_id=sponsor_id, cosponsor_ids=cosponsor_ids,
    cosponsors_count=cosponsors_count, withdrawn_cosponsors_count=withdrawn_cosponsors_count,
    committee_ids=committee_ids, related_bill_ids=related_bill_ids, enacted_as.congress=enacted_as.congress,
    enacted_as.number=enacted_as.number))

  # return_obj(as, query(url, args, ...))
  give_cg(as, url, "", args, ...)
}

ll <- function(x) {
  if (!is.null(x)) {
    if (!is.logical(x)) {
      stop(deparse(substitute(x)), " must be logical",
                             call. = FALSE)
    }
    if (x) {
      tolower(x)
    } else {
      x
    }
  }
}
