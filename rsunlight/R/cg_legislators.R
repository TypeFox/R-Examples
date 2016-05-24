#' Search for legislators.
#'
#' @template getleg
#'
#' @param latitude latitude of coordinate
#' @param longitude longitude of coordinate
#' @param zip zip code to search
#' @template cg
#' @template cg_query
#'
#' @return List of output fields.
#' @export
#' @examples \dontrun{
#' cg_legislators(last_name = 'Pelosi')
#' cg_legislators(party = 'D')
#' cg_legislators(facebook_id = 'mitchmcconnell')
#' cg_legislators(latitude = 35.778788, longitude = -78.787805)
#' cg_legislators(zip = 77006)
#'
#' # Output a list
#' cg_legislators(last_name = 'Pelosi', as='list')
#' # Output an httr response object, for debugging purposes
#' cg_legislators(last_name = 'Pelosi', as='response')
#'
#' # Pagination
#' cg_legislators(party = 'D', per_page=2)
#'
#' # Curl debugging
#' library('httr')
#' cg_legislators(party = 'D', config=verbose())
#' cg_legislators(party = 'D', config=timeout(0.1))
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_legislators(party = c('D', 'R'))
#' cg_legislators(last_name = c('Pelosi', 'Merkley'))
#' }

cg_legislators <- function(title=NULL, first_name=NULL, middle_name=NULL,
    last_name=NULL, name_suffix=NULL, nickname=NULL, party=NULL, state=NULL,
    state_name=NULL, state_rank=NULL, district=NULL, in_office=NULL, chamber=NULL,
    gender=NULL, phone=NULL, fax=NULL, office=NULL, website=NULL, contact_form=NULL,
    email=NULL, congress_office=NULL, bioguide_id=NULL, ocd_id=NULL, thomas_id=NULL,
    lis_id=NULL, crp_id=NULL, icpsr_id=NULL, votesmart_id=NULL, fec_ids=NULL,
    govtrack_id=NULL, congresspedia_url=NULL, twitter_id=NULL, youtube_id=NULL,
    facebook_id=NULL, senate_class=NULL, term_start=NULL, term_end=NULL, birthday=NULL,
    latitude = NULL, longitude = NULL, zip = NULL, query=NULL, fields=NULL, page=1, per_page=20,
    order=NULL, key = NULL, as = 'table', ...) {

  key <- check_key(key)
  if (!is.null(latitude) | !is.null(latitude) | !is.null(zip)) {
    url <- paste0(cgurl(), "/legislators/locate")
    if (!all(is.null(title),is.null(first_name),is.null(middle_name),is.null(last_name),
        is.null(name_suffix),is.null(nickname),is.null(party),is.null(state),is.null(state_name),
        is.null(state_rank),is.null(district),is.null(in_office),is.null(chamber),is.null(gender),
        is.null(phone),is.null(fax),is.null(office),is.null(website),is.null(contact_form),
        is.null(email),is.null(congress_office),is.null(bioguide_id),is.null(ocd_id),
        is.null(thomas_id),is.null(lis_id),is.null(crp_id),is.null(icpsr_id),is.null(votesmart_id),
        is.null(fec_ids),is.null(govtrack_id),is.null(congresspedia_url),is.null(twitter_id),
        is.null(youtube_id),is.null(facebook_id),is.null(senate_class),is.null(term_start),
        is.null(term_end),is.null(birthday)))
      stop("If latitude, longitude, or zip are used, all other parameters must be NULL", call. = FALSE)
    if (!is.null(latitude) & !is.null(latitude) & is.null(zip)) {
      stopifnot(is.null(zip))
      args <- sc(list(apikey=key,latitude=latitude,longitude=longitude,per_page=per_page,page=page,fields=fields,order=order))
    } else if (!is.null(zip)) {
      stopifnot(is.null(latitude),is.null(longitude))
      args <- sc(list(apikey=key,zip=zip,per_page=per_page,page=page,fields=fields,query=query,order=order))
    }
  } else {
    url <- paste0(cgurl(), "/legislators")
    args <- sc(list(apikey=key,title=title,first_name=first_name,middle_name=middle_name,
        last_name=last_name,name_suffix=name_suffix,nickname=nickname,party=party,state=state,
        state_name=state_name,state_rank=state_rank,district=district,in_office=in_office,
        chamber=chamber,gender=gender,phone=phone,fax=fax,office=office,website=website,
        contact_form=contact_form,email=email,congress_office=congress_office,
        bioguide_id=bioguide_id,ocd_id=ocd_id,thomas_id=thomas_id,lis_id=lis_id,crp_id=crp_id,
        icpsr_id=icpsr_id,votesmart_id=votesmart_id,fec_ids=fec_ids,govtrack_id=govtrack_id,
        congresspedia_url=congresspedia_url,twitter_id=twitter_id,youtube_id=youtube_id,
        facebook_id=facebook_id,senate_class=senate_class,term_start=term_start,term_end=term_end,
        birthday=birthday,per_page=per_page,page=page,fields=fields,query=query,order=order))
  }

  give_cg(as, url, "", args, ...)
}
