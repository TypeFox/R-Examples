#' YEARMODA to Date conversion
#' 
#' Converts dates in YEARMODA format to R date format
#' 
#' Converts YEARMODA to R date
#' 
#' @param YEARMODA Date in YEARMODA format (e.g. 20160206 for 6th Feb 2016)
#' @return Date object
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' 
#' YEARMODA2Date(20001205)
#' YEARMODA2Date(19901003)
#' 
#'  
#' @export YEARMODA2Date
YEARMODA2Date<-function(YEARMODA)
{if(is.numeric(YEARMODA))
  {Year<-trunc(YEARMODA/10000)
Month<-trunc((YEARMODA-Year*10000)/100)
Day<-YEARMODA-Year*10000-Month*100
return(ISOdate(Year,Month,Day))} 
}
