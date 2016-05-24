#' Not sure how this differs from their other API...
#'
#' @export
#' @param genre Book genre
#' @param title Book title
#' @param aufirst First author
#' @param aulast Last author
#' @param date Date of publication
#' @param spage Start page
#' @param issue Issue number
#' @param version One of 0.1 or 1.0
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_openurl(
#' 	genre="book",
#' 	title="Manual+of+North+American+Diptera",
#' 	aufirst="Samuel Wendell",
#' 	aulast="Williston",
#' 	date=1908,
#' 	spage=16)
#'
#' bhl_openurl(genre="book", title="Manual+of+North+American+Diptera",
#'    aufirst="Samuel Wendell", aulast="Williston", date=1908, spage=16)
#'
#' bhl_openurl(genre="book", title="Manual+of+North+American+Diptera",
#'    aufirst="Samuel Wendell", aulast="Williston", date=1908, spage=16, as='xml')
#' }

bhl_openurl <- function(genre = NULL, title = NULL, aufirst = NULL, aulast = NULL,
	date = NULL, spage = NULL, issue = NULL, version = 0.1, as = "list",
  key = NULL, ...)
{

  if (version == "1.0") {
    url_ver <- "z39.88-2004"
  } else {
    url_ver <- NULL
  }
	args <- bhlc(list(genre = genre, title = title, aufirst = aufirst,
							 date = date, spage = spage, issue = issue, url_ver = url_ver,
							 apikey = check_key(key), format = as_f(as)))
	if (length(args) == 0) args <- NULL
	out <- GET("http://www.biodiversitylibrary.org/openurl", query = args, ...)
	stop_for_status(out)
	tt <- content(out, as = "text")
  switch(as, json = tt, xml = tt, list = fjson(tt), table = todf(tt))
}
