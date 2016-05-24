#' List the unique names.
#'
#' By using the startrow and batchsize parameters appropriately, you can
#'    pull the list all at once, or in batches (i.e. 1000 names at a time).
#'    Names both with and without NameBank identifiers are returned.
#'
#' @export
#'
#' @param startrow first name to return (if using as an offset)
#' @param batchsize number of names to return (numeric)
#' @param startdate (optional) start date of range between which to count names (date)
#' @param enddate (optional) end date of range between which to count names (date)
#' @inheritParams bhl_getcollections
#'
#' @examples \dontrun{
#' bhl_namelist(startrow=1, batchsize=99, startdate='10/15/2009', enddate='10/16/2009')
#' bhl_namelist(startrow=1, batchsize=5, startdate='10/15/2009', enddate='10/31/2009',
#'    as='json')
#' }

bhl_namelist <- function(startrow = NULL, batchsize = NULL, startdate = NULL,
  enddate = NULL, as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "NameList", apikey = check_key(key), format = as_f(as), startrow = startrow,
                       batchsize = batchsize, startdate = startdate, enddate = enddate))
  bhl_GET(as, args, ...)
}
