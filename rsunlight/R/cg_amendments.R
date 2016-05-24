#' Search for data on amendments.
#'
#' @export
#'
#' @param amendment_id The unique ID for this amendment. Formed from the amendment_type, number,
#' and congress.
#' @param amendment_type The type for this amendment. For the amendment 'H.Amdt. 10', the
#' amendment_type represents the 'H.Amdt.' part. Amendment types can be either hamdt or samdt.
#' @param number The number for this amendment. For the amendment 'H.Amdt. 10', the number is 10.
#' @param congress The Congress in which this amendment was introduced. For example, amendments
#' introduced in the '113th Congress' have a congress of 113.
#' @param chamber The chamber in which the amendment was introduced.
#' @param house_number If the amendment was introduced in the House, this is a relative amendment
#' number, scoped to the bill or treaty the House it relates to. How this number gets assigned is
#' complicated and involves multiple institutions within the House and the Library of Congress.
#' You can read the gory details if you want, but this number will usually do the job of connecting
#' to data from the House Clerk's Electronic Voting System.
#' @param introduced_on The date this amendment was introduced.
#' @param last_action_at The date or time of the most recent official action on the amendment.
#' Often, there are no official actions, in which case this field will be set to the value of
#' introduced_on.
#' @param amends_bill_id If this amendment relates to a bill, this field is the ID of the
#' related bill.
#' @param amends_treaty_id If this amendment relates to a treaty, this field is the ID of the
#' related treaty. Treaty IDs are of the form treatyX-Y, where X is the treaty's number, and Y
#' is the Congress the treaty is being considered in.
#' @param amends_amendment_id If this amendment amends an amendment, this field is the ID of the
#' amended amendment.
#' @param sponsor_type Whether the amendment is sponsored by a 'person' or a 'committee'.
#' @param sponsor_id If the sponsor_type is 'person', this will be that legislator's bioguide ID.
#' If the sponsor_type is 'committee', this will be that committee's ID.
#'
#' @template cg
#' @template cg_query
#' @examples \dontrun{
#' cg_amendments()
#' cg_amendments(chamber='house', congress=113)
#' cg_amendments(sponsor_type='committee', sponsor_id='HSRU')
#' cg_amendments(amends_bill_id='hr624-113')
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_amendments(chamber = c('house', 'senate'), per_page=2)
#' }

cg_amendments <- function(amendment_id=NULL, amendment_type=NULL, number=NULL, congress=NULL,
  chamber=NULL, house_number=NULL, introduced_on=NULL, last_action_at=NULL, amends_bill_id=NULL,
  amends_treaty_id=NULL, amends_amendment_id=NULL, sponsor_type=NULL, sponsor_id=NULL,
  query=NULL, fields=NULL, page=1, per_page=20, order=NULL,
  key = NULL, as = 'table', ...) {

  key <- check_key(key)
  url <- 'https://congress.api.sunlightfoundation.com/amendments'
  args <- sc(list(apikey=key, amendment_id=amendment_id, amendment_type=amendment_type,
                          number=number, congress=congress, chamber=chamber, house_number=house_number,
                          introduced_on=introduced_on, last_action_at=last_action_at,
                          amends_bill_id=amends_bill_id,amends_treaty_id=amends_treaty_id,
                          amends_amendment_id=amends_amendment_id, sponsor_type=sponsor_type,
                          sponsor_id=sponsor_id, query=query, per_page=per_page, page=page,
                          fields=fields, order=order))
  give_cg(as, url, "", args, ...)
}
