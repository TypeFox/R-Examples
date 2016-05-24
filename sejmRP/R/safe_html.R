#' Safe html scrapping
#'
#' Function \code{safe_html} tries to download the URL several times.
#' 
#' @details
#' Function \code{safe_html} performes 10 (by default) attempts to download the URL
#' and waits 60sec (by default) after each failure
#'
#' @usage safe_html(page, time = 60, attempts = 10)
#'
#' @param page requested URL
#' @param time sleep interval after each failure
#' @param attempts max number of tries (if there is a problem with connection)
#'
#' @return character vector
#'
#' @examples
#' \dontrun{
#' page <- paste0('http://www.sejm.gov.pl/Sejm7.nsf/',
#'                'wypowiedz.xsp?posiedzenie=15&dzien=1&wyp=008')
#' safe_html(page)}
#'
#' @author Przemyslaw Biecek
#'
#' @export
#' @importFrom xml2 read_html
#'

safe_html <- function(page, time = 60, attempts = 10) {
  stopifnot(is.character(page))
  stopifnot(is.numeric(time))
  stopifnot(is.numeric(attempts))
  
  repeat({
    attempts <- attempts - 1
    pageH <- try(read_html(page), silent = TRUE)
    if (class(pageH)[1] != "try-error")
      break()
    cat("No connection, trying again: ", page, "\n")
    Sys.sleep(time)
    if (attempts < 0) {
      stop("No internet connection")
    }
  })

  return(pageH)
} 
