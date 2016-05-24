#' Plot Estimations of Storm Activity Preindex
#'
#'This function is used to visualize the preindex estimation for each station.
#'Users can specify the number of graphs per page, the default is set to 2 graphs per page.
#'
#'@param x estimation of preindex from \code{SAIndex}
#'@param Title title of preindexplot
#'@param start start date of records for magnetic activities. See examples
#'@param end end date of records for magnetic activities. See examples
#'@param n.station  number of stations
#'@param graphs.per.page how many graphs combined in one plot page. Default number is 2
#'@param station.names NULL (indicating no station names) or a vector strings for station names
#'@param ... additional arguments.
#'
#'@details This function is used to visualize the preindex. The function plots the preindex estimation from \code{SAIndex}. See example below. 
#'@author Inga Maslova, Wen Xiao
#'@seealso \code{\link{SAIndex}}
#'
#'@examples
#'\dontrun{
#'## data sets.
#'data(record)
#'## Colatitude and longitude of geomagnetic observatories.
#'coord=matrix(c(124.43, 19.23, 53.77,140.18,68.68,202.00,71.89,293.85),nrow=2,ncol=4)
#'
#'
#'## estimation of one iWISA and one preindex for each stations.
#'data<- SAIndex(record, coord, wf="la8")
#'
#'## generate datetime for one week period
#'start.date="2001-3-1"
#'end.date="2001-4-30"
#'
#'preindexplot(data, Title="Preindex of stations", start = start.date,end=end.date,
#'n.station=4, graphs.per.page=2, station.names=c("HER","KAK","HON","SJG"))
#'}
#'@keywords Preindexplot
#'@export
#'@importFrom graphics axis.POSIXct par plot 

preindexplot <- function(x, Title=NULL, start = NULL,end=NULL,n.station=NULL,
                          graphs.per.page=2, station.names=NULL,...){
  N<-nrow(data.frame(x))
  indexlist<-matrix(0, ncol=n.station, nrow=N)
  indexlist<-data.frame(x[-1])
  colnames(indexlist) <-station.names
  # paste("Station", station.names, sep = " ")

  days<-seq(as.Date(start), as.Date(end), by = "day")

  dates.value<- data.frame(Date=as.Date(rep(days,each=1440)),
                           Time=format(as.difftime(0:1439, units = "mins") +
                                         as.POSIXct(start), "%H:%M"))
  index.value<-cbind(dates.value, indexlist)
  index.value$datetime <- as.POSIXct(paste(index.value$Date,index.value$Time),tz="UCT")
  index.final<-index.value[,-(1:2)]
  ### Making a series of plots that proceed by a click
  #build the data frame
  mr <- par()$mar
  Ngraphs <- min(n.station, graphs.per.page)
  par(ask=TRUE)
  par(mfrow=c(Ngraphs, 1))
  col=1:n.station
  n=nrow(index.final)
  for(i in 1:n.station){
    plot(index.final$datetime, index.final[,i],xaxt = "n", type="l", col=col[i], main=Title,
         xlab="DateTime", ylab=paste("Station",colnames(index.final[i])))
    r <- as.POSIXct(round(range(index.final$datetime), "hours"))
    axis.POSIXct(1, at=seq(r[1], r[2], by = "days"),  format="%m/%d"
    )
  }
}
