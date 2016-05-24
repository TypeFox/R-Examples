#' Get the last date from a phenology record
#' 
#' When looking at multi-year phenology records, it is normally obvious in
#' which year bloom occurred last. Determining this with an automated
#' procedure, however, is a bit tricky, when the range of phenological dates
#' spans across a calendar year transition. This function finds the latest
#' phenological date of the record. This is the date before the longest
#' phenological date gap.
#' 
#' 
#' @param dates numeric vector of Julian dates (days of the year)
#' @param first boolean variable that can be set to TRUE to get the first, not
#' the last, date of the phenology record.
#' @return the latest (earliest) date of the series, under the assumption that
#' the longest period without bloom can be interpreted as separating the
#' phenological seasons. This should be a reasonable assumption in most cases.
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' 
#' get_last_date(c(1,3,6,8,10,25))
#' get_last_date(c(345,356,360,365,2,5,7,10))
#' get_last_date(c(345,356,360,365,2,5,7,10),first=TRUE)
#' 
#'  
#' @export get_last_date
get_last_date<-function(dates,first=FALSE)
{pdates<-dates[order(dates)]
gaps<-c(pdates[2:length(pdates)]-pdates[1:(length(pdates)-1)],
        which(pdates[1]==c(180:366,1:179))-which(pdates[length(pdates)]==c(180:366,1:179)))
gaps[which(gaps<0)]<-365+gaps[which(gaps<0)]
last<-max(pdates[which(gaps==max(gaps))])
if(which(gaps==max(gaps))==length(pdates)) firsty<-pdates[1] else firsty<-pdates[which(gaps==max(gaps))+1]
#first<-max(pdates[which(gaps==max(gaps))])
if(first) return(firsty) else return(last)
}
