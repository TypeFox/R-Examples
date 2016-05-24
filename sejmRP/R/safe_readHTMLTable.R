#' Safe html table scrapping
#'
#' Function \code{safe_readHTMLTable} tries to download the table from given URL several times.
#' 
#' @details
#' Function \code{safe_readHTMLTable} performes 10 (by default) attempts to download the URL
#' and waits 60sec (by default) after each failure
#'
#' @usage safe_readHTMLTable(..., time = 60, attempts = 10)
#'
#' @param ... arguments that will be passed to readHTMLTable
#' @param time sleep interval after each failure
#' @param attempts max number of tries (if there is a problem with connection)
#'
#' @return character vector
#'
#' @examples
#' \dontrun{
#' page <- paste0('http://www.sejm.gov.pl/Sejm7.nsf/',
#'                'posiedzenie.xsp?posiedzenie=99&dzien=2')
#' safe_readHTMLTable(page)}
#'
#' @author Przemyslaw Biecek
#'
#' @export
#' @importFrom XML readHTMLTable
#'

safe_readHTMLTable <- function(..., time = 60, attempts = 10) {
  stopifnot(is.numeric(time))
  stopifnot(is.numeric(attempts))
  
  repeat({
    attempts <- attempts - 1
    pageH <- try(readHTMLTable(...), silent = TRUE)
    if (class(pageH)[1] != "try-error")
      break()
    cat("No connection, trying again: \n")
    Sys.sleep(time)
    if (attempts < 0) {
      stop("No internet connection")
    }
  })

  return(pageH)
} 
