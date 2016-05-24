#' Days since 14. Oct 1582
#'
#' Count the number of days since October 14, 1582.
#'
#' @usage xpssJDate()
#' @return Returns days since 14 Oct 1582
#' @author Bastian Wiessner
#' @examples
#' xpssJDate()
#' @export

xpssJDate <- function(){
  
  gregorian <- as.Date(x="1582-10-14")
  x <- difftime(time2=Sys.Date(),time1=gregorian)
  x <- as.numeric(paste(x*-1))
  
  return(x)
}

