#' Search for titles and items in BHL.
#'
#' Search criteria includes title, author last name, volume, edition, year of
#'    publication, subject, language code, and collection identifier. Valid
#'    language codes and collection identifiers can be obtained from the
#'    getlanguages and getcollections functions. If year of publication is
#'    specified, it should be a 4-digit year. To execute a search, you must
#'    supply at least a title, author last name, or collection identifier.
#'
#' @export
#'
#' @param title string to search for in the title (character)
#' @param lname last name to search for (character)
#' @param volume volume to search for (numeric)
#' @param edition edition to search for (character)
#' @param year year to search for, four characters, e.g, 1970 (numeric)
#' @param collectionid collection identifier to search for (numeric)
#' @param language language to search for (character)
#' @inheritParams bhl_getcollections
#'
#' @note Use \code{\link{bhl_getcollections}} or \code{\link{bhl_getlanguages}} to get
#' acceptable terms
#' @examples \dontrun{
#' bhl_booksearch(title='Selborne', lname='White', volume=2, edition='new', year=1825,
#'    collectionid=4, language='eng')
#' bhl_booksearch(title='evolution', year=2000, as='json')
#' bhl_booksearch('evolution', year=2000, as='xml')
#' bhl_booksearch('evolution', year=2000, as="list")
#' }
bhl_booksearch <- function(title = NULL, lname = NULL, volume = NULL,
  edition = NULL, year = NULL, collectionid = NULL, language = NULL,
  as='table', key = NULL, ...) {

  args <- bhlc(list(op = "BookSearch", apikey = check_key(key), format = as_f(as),
                       title=title, lname=lname, volume=volume,
                       edition=edition, year=year, collectionid=collectionid,
                       language=language))
  bhl_GET(as, args, ...)
}
