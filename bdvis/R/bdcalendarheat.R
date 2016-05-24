#' Calendar heat map of biodiversity data
#' 
#' Produces a heat map {\url{https://en.wikipedia.org/wiki/Heat_map}}
#' representing the distribution of records in time.
#' 
#' The calendar heat map is a matrix-like plot where each cell represents a 
#' unique date, and the color the cell is painted with shows the amount of 
#' records that have that particular date. Rows are weekdays and columns are 
#' week numbers, each year having its own "panel".
#' 
#' @import sqldf
#' @importFrom stats na.omit
#' @param indf input data frame containing biodiversity data set
#' @param title title custome title for the plot
#' @examples \dontrun{
#' bdcalendarheat(inat)
#' }
#' @family Temporal visualizations
#' @export
bdcalendarheat <- function(indf=NA,title=NA){
  if(is.na(title)){
    title="number of records"
  }
  indf$Date_collected = as.Date(indf$Date_collected,"%Y-%m-%d")
  dat=sqldf("select Date_collected, count(*) as recs from indf group by Date_collected")
  #dat=dat[2:dim(dat)[1],]
  dat=na.omit(dat)
  Year = as.numeric(strftime(as.Date(dat$Date_collected,na.rm=T), format = "%Y"))
  CurrentYear = as.numeric(strftime(as.Date(Sys.Date()), format = "%Y"))
  if(max(Year)>CurrentYear){
    dat=dat[which(Year <= CurrentYear ),]
  }
  Year = as.numeric(strftime(as.Date(dat$Date_collected,na.rm=T), format = "%Y"))
  if(max(Year)-min(Year) > 6) {
    dat=dat[which(Year > (max(Year)-6) ),]
  }
  calendarHeat(dat$Date_collected, dat$recs, varname=title)
}
