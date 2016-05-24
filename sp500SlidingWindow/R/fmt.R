#' Format a Number with Commas and No Decimals
#'
#' @author George Fisher
#'
#' #@description
#'
#' #@details
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dat = 123456.654
#' fmt(dat)
#' }
#'
#' @return formatted string
#'
#' @param dat number to be formatted
#'
fmt <- function(dat) {
    return(format(dat, big.mark = ",", scientific = FALSE, digits = 0))
}
