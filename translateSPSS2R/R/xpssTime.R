#' Minutes since 14. Oct 1582
#'
#'  R implementation of the SPSS \code{$TIME} function. Count of the number of minutes since October 14, 1582 (the first day of the Gregorian calendar).
#'
#' @usage xpssTime()
#' @return Returns minutes since 14 Oct 1582
#' @author Bastian Wiessner
#' @export

xpssTime <- function(){
  
  gregorian <- as.Date(x="1582-10-14")
  x <- difftime(time2=Sys.Date(),time1=gregorian,units="sec")
  x <- as.numeric(paste(x*-1))
  #manuelle anpassung
  x <- x+56950
  return(x)
}
