#' Return the number of unique names found on pages in BHL.
#'
#' Names both with and without NameBank identifiers are counted.
#'
#' BEWARE: this API call can take a long time. They are likely working on speeding up the
#' service, but slow for now.
#'
#' @export
#' @param startdate start date of range between which to count names (optional)
#' @param enddate end date of range between which to count names (optional)
#' @inheritParams bhl_getcollections
#'
#' @examples \dontrun{
#' bhl_namecount(startdate = '12/25/2009', enddate = '12/27/2009')
#' bhl_namecount(startdate = '10/15/2009', enddate = '10/17/2009', as='json')
#' bhl_namecount(startdate = '10/15/2009', enddate = '10/17/2009', as='xml')
#' }

bhl_namecount <- function(startdate = NULL, enddate = NULL, as='table', key = NULL, ...)
{
  args <- bhlc(list(op = "NameCount", apikey = check_key(key), format = as_f(as),
                       startdate=startdate, enddate=enddate))
  bhl_GET(as, args, ...)
}
