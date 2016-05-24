#' Search for data on votes.
#'
#' @export
#'
#' @param roll_id A unique identifier for a roll call vote. Made from the first letter of the
#' chamber, the vote number, and the legislative year.
#' @param chamber The chamber the vote was taken in. 'house' or 'senate'
#' @param number The number that vote was assigned. Numbers reset every legislative year.
#' @param year The 'legislative year' of the vote. This is not quite the same as the calendar year
#' - the legislative year changes at noon EST on January 3rd. A vote taken on January 1, 2013 has
#' a 'legislative year' of 2012.
#' @param congress The Congress this vote was taken in.
#' @param voted_at The time the vote was taken.
#' @param vote_type The type of vote being taken. This classification is imperfect and unofficial,
#' and may change as we improve our detection. Valid types are 'passage', 'cloture', 'nomination',
#' 'impeachment', 'treaty', 'recommit', 'quorum', 'leadership', and 'other'.
#' @param roll_type The official description of the type of vote being taken.
#' @param required The required ratio of Aye votes necessary to pass the legislation. A value of
#' '1/2' actually means more than 1/2. Ties are not possible in the Senate (the Vice President
#' casts a tie-breaker vote), and in the House, a tie vote means the vote does not pass.
#' @param result The official result of the vote. This is not completely standardized (both
#' 'Passed' and 'Bill Passed' may appear). In the case of a vote for Speaker of the House, the
#' result field contains the name of the victor.
#' @param bill_id If a vote is related to a bill, the bill's ID.
#' @param nomination_id If a vote is related to a nomination, the nomination's ID.
#' @param callopts Curl options to be passed on to httr::GET
#' @param ... See Details. You can pass on additional parameters.
#' @param fields You can request specific fields by supplying a vector of fields names. Many fields
#' are not returned unless requested. If you don't supply a fields parameter, you will get the
#' most commonly used subset of fields only. To save on bandwidth, parsing time, and confusion,
#' it's recommended to always specify which fields you will be using.
#' @param key your SunlightLabs API key; loads from .Rprofile
#' @param as (character) One of table (default), list, or response (httr response object)
#' @param page Page to return. Default: 1. You can use this in combination with the
#' per_page parameter to get more than the default or max number of results per page.
#' @param per_page Number of records to return. Default: 20. Max: 50.
#' @param order Sort results by one or more fields with the order parameter. order is
#' optional, but if no order is provided, the order of results is not guaranteed to be predictable.
#' Append \code{__asc} or \code{__desc} to the field names to control sort direction. The default
#' direction is \code{desc}, because it is expected most queries will sort by a date. Any field
#' which can be used for filtering may be used for sorting. On full-text search endpoints (URLs
#' ending in \code{/search}), you may sort by score to order by relevancy.
#'
#' @template cg_query
#'
#' @details Two parameters can be passed on that vary with vote and/or party plus vote:
#' \itemize{
#'  \item breakdown.total.[vote] The number of members who cast [vote], where [vote] is a valid
#'  vote as defined above.
#'  \item breakdown.party.[party].[vote] The number of members of [party] who cast [vote],
#'  where [party] is one of 'D', 'R', or 'I', and [vote] is a valid vote as defined above.
#' }
#'
#' @examples \dontrun{
#' cg_votes(chamber='senate', order='voted_at')
#' cg_votes(query='guns')
#' cg_votes(voter_ids.A000055__exists=TRUE)
#'
#' # Pass in more than one value for a parameter
#' cg_votes(chamber = c('house', 'senate'))
#' }

cg_votes <- function(roll_id=NULL, chamber=NULL, number=NULL, year=NULL, congress=NULL,
  voted_at=NULL, vote_type=NULL, roll_type=NULL, required=NULL, result=NULL, bill_id=NULL,
  nomination_id=NULL, query=NULL, fields=NULL, page=1, per_page=20, order=NULL,
  key = NULL, as = 'table', callopts = list(), ...) {

  key <- check_key(key)
  args <- sc(list(apikey=key,roll_id=roll_id, chamber=chamber, number=number, year=year,
      congress=congress, voted_at=voted_at, vote_type=vote_type, roll_type=roll_type,
      required=required, result=result, bill_id=bill_id, nomination_id=nomination_id,
      per_page=per_page, page=page, fields=fields, order=order, ...))

  # return_obj(as, query(paste0(cgurl(), "/votes"), args, callopts, ...))
  give_cg(as, cgurl(), "/votes", args, callopts)
}
