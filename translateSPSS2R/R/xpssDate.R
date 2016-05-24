#' Today's date 
#'
#' R implementation of the SPSS \code{$date} system variable. 
#' 
#' xpssDate provides the SPSS format of dates.
#' 
#' @usage xpssDate()
#' @return Current date as "dd-mmm-yy". Day and year are numeric values with the length two, month is a character string with the length three.
#' @author Bastian Wiessner
#' @examples 
#' xpssDate()
#' @export

xpssDate <- function(){

  (today <- Sys.Date())
  today <- toupper(today)
  
  return(today)
}

as.POSIXct
