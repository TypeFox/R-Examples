#' Search for parts of books in BHL, such as articles, chapters, or treatments.
#' Search criteria includes title, container (journal or book title), author,
#' date of publication, volume, series, and issue.
#'
#' To execute a search, you must supply at least a title or author.
#'
#' The metadata returned by this method includes Part Identifier, Part URL, Item ID,
#' Page ID for the start page, Genre, Title, Container Title, Publication Details,
#' Volume, Series, Issue, Date, Page Range, Language, rights information, authors,
#' keywords, identifiers, pages, and related parts. For more information, see the
#' "Data Elements" section of this documentation.
#'
#' @export
#' @param title Title of the work
#' @param containerTitle Container title of the work
#' @param author Author of the work
#' @param date Date of the work
#' @param volume Volume of the work
#' @param series Series of the work
#' @param issue Issue of the work
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_partsearch(title='Critical approach to the definition of Darwinian units')
#' bhl_partsearch(author='Charles Darwin')
#' }

bhl_partsearch <- function(title=NULL, containerTitle=NULL, author=NULL, date=NULL,
  volume=NULL, series=NULL, issue=NULL, as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "PartSearch", apikey = check_key(key), format = as_f(as),
                       title = title, containerTitle=containerTitle, author=author,
                       date=date, volume=volume, series=series, issue=issue))
  bhl_GET(as, args, ...)
}
