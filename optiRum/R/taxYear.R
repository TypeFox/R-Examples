#' Returns the UK financial tax year for a given date
#'
#' Base don UK tax year April 6 - April 5, this returns the year (YYYY) 
#' the tax period covers. Tax Year start date can be overriden.
#'
#' @param date    Date to be checked
#' @param start   Provide the month & day that will be used as the first tax day (mm-dd)
#' 
#' @return year   The financial year
#'
#' @keywords financial tax
#' @family tax
#' @export
#' 
#' @examples
#' # single set of values
#' taxYear(Sys.Date()) 
#' 
#' # vector of values
#' taxYear(seq(Sys.Date(),by=1,length=500))
#' 

taxYear<-function(date=Sys.Date(),start="04-06"){
  
  calendaryear<-year(date)
  
  currenttaxboundary<-as.Date(paste(calendaryear,start,sep="-"))
  
  afterstart<-date>=currenttaxboundary
  
  year<-calendaryear-1L+afterstart
    
  return(year)
}

