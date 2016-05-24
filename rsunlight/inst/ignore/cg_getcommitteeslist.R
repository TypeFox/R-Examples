#' Get list of all committees for a given chamber along with their subcommittees.
#'
#' @import httr
#' @param chamber House, Senate, or Joint
#' @template cg
#' @return List of all committees in the specified chamber with their
#'    subcommittees (but not memberships due to size of response).
#' @export
#' @examples \dontrun{
#' cg_getcommitteeslist(chamber = 'Joint')
#' }
cg_getcommitteeslist <- function(chamber = NULL,
  key=getOption("SunlightLabsKey", stop("need an API key for Sunlight Labs")),
  callopts = list())
{
  url = "http://services.sunlightlabs.com/api/committees.getList.json"
  args <- suncompact(list(apikey = key, chamber = chamber))
  tt <- GET(url, query=args, callopts)
  stop_for_status(tt)
  out <- content(tt, as = "text")
  res <- fromJSON(out, simplifyVector = FALSE)$response
  class(res) <- "cg_getcommitteelist"
  return( res )
}
