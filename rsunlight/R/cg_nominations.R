#' Search for data on nominations
#'
#' @export
#'
#' @param nomination_id The unique identifier for this nomination, taken from the Library of
#' Congress. Of the form 'PN[number]-[congress]'.
#' @param congress The Congress in which this nomination was presented.
#' @param number The number of this nomination, taken from the Library of Congress. Can
#' occasionally contain hyphens, e.g. 'PN64-01'.
#' @param received_on The date this nomination was received in the Senate.
#' @param last_action_at The date this nomination last received action. If there are no official
#' actions, then this field will fall back to the value of received_on.
#' @param organization The organization the nominee would be appointed to, if confirmed.
#' @param committee_ids An array of IDs of committees that the nomination has been referred to
#' for consideration.
#' @param nominees An array of objects with fields (described below) about each nominee.
#' Nominations for civil posts tend to have only one nominee. Nominations for military posts tend
#' to have batches of multiple nominees. In either case, the nominees field will always be an array.
#' @param nominees.position The position the nominee is being nominated for.
#' @param nominees.state The which state in the United States this nominee hails from. This field
#' is only available for some nominees, and never for batches of multiple nominees.
#'
#' @template cg
#' @template cg_query
#' @examples \dontrun{
#' cg_nominations(order='received_on')
#' cg_nominations(committee_ids='SSAS')
#' cg_nominations(organization='Privacy and Civil Liberties Oversight Board')
#' cg_nominations(query='Petraeus')
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_nominations(party = c('PN1873-111', 'PN604-112'))
#' }

cg_nominations <- function(nomination_id=NULL, congress=NULL, number=NULL, received_on=NULL, last_action_at=NULL,
  organization=NULL, committee_ids=NULL, nominees=NULL, nominees.position=NULL, nominees.state=NULL,
  query=NULL, fields=NULL, page=1, per_page=20, order=NULL,
  key = NULL, as = 'table', ...) {

  key <- check_key(key)
  url <- 'https://congress.api.sunlightfoundation.com/nominations'
  args <- sc(list(apikey=key,nomination_id=nomination_id, congress=congress,
    number=number, received_on=received_on, last_action_at=last_action_at, organization=organization,
    committee_ids=committee_ids, nominees=nominees, nominees.position=nominees.position,
    nominees.state=nominees.state, query=query, per_page=per_page, page=page, fields=fields,
    order=order))
  give_cg(as, url, "", args, ...)
}
