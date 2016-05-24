#' Search for data on hearings
#'
#' @export
#'
#' @param committee_id (numeric) The ID of the committee holding the hearing.
#' @param occurs_at (character) The time the hearing will occur. format: "%Y-%m-%dT%H:%M:%SZ"
#' e.g. "2014-05-01T13:00:00Z"
#' @param congress (numeric) The number of the Congress the committee hearing is taking
#' place during.
#' @param chamber (character) The chamber ('house', 'senate', or 'joint') of the committee
#' holding the hearing.
#' @param dc (logical) Whether the committee hearing is held in DC \code{TRUE} or in the
#' field \code{FALSE}
#' @param bill_ids (numeric) The IDs of any bills mentioned by or associated with the hearing.
#' @param hearing_type (character) (House only) The type of hearing this is. Can be: 'Hearing',
#' 'Markup', 'Business Meeting', 'Field Hearing'.
#'
#' @template cg
#' @template cg_query
#'
#' @return url (character) (House only) A permalink to that hearing's description on that
#' committee's official website.
#' @return description (character) A description of the hearing.
#' @return room (character) If the hearing is in DC, the building and room number the hearing
#' is in. If the hearing is in the field, the address of the hearing.
#' @return hearing_id (character) (House only) A permalink to that hearing's description on that
#' committee's official website.
#'
#' @references
#' \url{https://sunlightlabs.github.io/congress/hearings.html}
#'
#' @examples \dontrun{
#' cg_hearings(chamber='house', dc=TRUE)
#' cg_hearings(query='children')
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_hearings(chamber = c('house', 'senate'))
#' }

cg_hearings <- function(committee_id=NULL, occurs_at=NULL, congress=NULL, chamber=NULL,
  dc=NULL, bill_ids=NULL, hearing_type=NULL, query=NULL, fields=NULL, page=1, per_page=20, order=NULL,
  key = NULL, as = 'table', ...) {

  key <- check_key(key)
  args <- sc(list(apikey=key, committee_id=committee_id, occurs_at=occurs_at,
      congress=congress, chamber=chamber, dc=getdc(dc), bill_ids=bill_ids,
      hearing_type=hearing_type, query=query, per_page=per_page, page=page, fields=fields, order=order))
  give_cg(as, cgurl(), "/hearings", args, ...)
}

getdc <- function(x){
  if (is.null(x)) {
    NULL
  } else if (x) {
    'true'
  } else  {
    'false'
  }
}
