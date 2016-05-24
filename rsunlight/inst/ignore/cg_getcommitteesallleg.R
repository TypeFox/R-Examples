#' Gets a list of all committees that a member serves on, including subcommittes.
#' 
#' @import httr
#' @template cg
#' @param bioguide_id legislator's bioguide_id
#' @return Committee details for all committees that the given member serves on.
#' @export
#' @examples \dontrun{
#' cg_getcommitteesallleg(bioguide_id = 'S000148')
#' }
cg_getcommitteesallleg <- function(bioguide_id = NULL,
    key=getOption("SunlightLabsKey", stop("need an API key for Sunlight Labs")),
    callopts = list()) 
{
  url = "http://services.sunlightlabs.com/api/committees.allForLegislator.json"
  args <- suncompact(list(apikey = key, bioguide_id = bioguide_id))
  tt <- GET(url, query=args, callopts)
  stop_for_status(tt)
  out <- content(tt, as = "text")
  res <- fromJSON(out, simplifyVector = FALSE)$response
  class(res) <- "cg_getcommitteesallleg"
  return( res )
}