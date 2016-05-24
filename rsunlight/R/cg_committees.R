#' @title Gets details (subcommittees + membership) for a committee by id.
#'
#' @description Names, IDs, contact info, and memberships of committees and subcommittees
#' in the House and Senate. All committee information is sourced from bulk data at
#' github.com/unitedstates, which in turn comes from official House and Senate sources.
#' Feel free to open a ticket with any bugs or suggestions. We only provide information
#' on current committees and memberships. For historic data on committee names, IDs,
#' and contact info, refer to the bulk data.
#'
#' @export
#' @param member_ids An array of bioguide IDs of legislators that are assigned to
#' this committee.
#' @param committee_id Official ID of the committee, as it appears in various official
#' sources (Senate, House, and Library of Congress).
#' @param chamber The chamber this committee is part of. 'house', 'senate', or 'joint'.
#' @param subcommittee Whether or not the committee is a subcommittee.
#' @param parent_committee_id If the committee is a subcommittee, the ID of its
#' parent committee.
#' @template cg
#' @template cg_query
#' @return Committee details including subcommittees and all members.
#' @examples \dontrun{
#' cg_committees(member_ids='L000551')
#' cg_committees(committee_id='SSAP')
#' cg_committees(committee_id='SSAP', fields='members')
#' cg_committees(chamber='joint', subcommittee=FALSE)
#' cg_committees(parent_committee_id='HSWM')
#'
#' # Disable pagination
#' cg_committees(per_page='all')
#'
#' # Output a list
#' cg_committees(member_ids='L000551', as='list')
#' # Output an httr response object, for debugging purposes
#' cg_committees(member_ids='L000551', as='response')
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_committees(committee_id = c('SSAP', 'SSGA01'))
#' }

cg_committees <-  function(member_ids = NULL, committee_id = NULL, chamber = NULL, subcommittee = NULL,
  parent_committee_id = NULL, query=NULL, fields = NULL, page = 1, per_page = 20, order = NULL,
  key = NULL, as = 'table', ...) {

  key <- check_key(key)
  if (!is.null(subcommittee)) {
    subcommittee <- ifelse(subcommittee, 'true', 'false')
  }
  fields <- paste0(fields, collapse = ",")
  args <- sc(list(apikey = key, member_ids = member_ids, committee_id = committee_id,
              chamber = chamber, subcommittee = subcommittee, fields = fields,
              parent_committee_id = parent_committee_id, page = page, per_page = per_page,
              query = query, order = order))
  give_cg(as, cgurl(), "/committees", args, ...)
}
