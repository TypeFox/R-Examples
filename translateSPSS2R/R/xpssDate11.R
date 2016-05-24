#' Todays date 
#' 
#'  R implementation of the SPSS \code{$date11} system variable. 
#'
#' Provides the SPSS format of dates
#'
#' @usage xpssDate11()
#' @return Current date as "dd-mmm-yyyy". Day is a numeric with the length two, month is a character string with the length three and year is a numeric value with the length four. 
#' @author Bastia Wiessner
#' @examples
#' xpssDate11()
#' @export
xpssDate11 <- function(){
  
  (today <- Sys.Date())
  today <- toupper(today)
  
  return(today)
}
