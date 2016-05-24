#' Return a list of titles associated with a subject.
#'
#' @export
#' @param subject the full or partial subject for which to search (character)
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getsubjecttitles('diptera')
#' bhl_getsubjecttitles('diptera', 'xml')
#' bhl_getsubjecttitles('diptera', 'json')
#' }

bhl_getsubjecttitles <- function(subject = NULL, as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "GetSubjectTitles", apikey = check_key(key), format = as_f(as),
                       subject = subject))
  bhl_GET(as, args, ...)
}
