#'Plot Global Storm Activity Index WISA
#'
#'Plot WISA, wavelet index of storm activity, by using ggplot2
#'
#'@param x estimation of WISA from \code{SAIndex}
#'@param type line width (default=1)
#'@param main, xlab, ylab graphical arguments, see \code{\link{par}}
#'@param start start date of records for magnetic activities.See examples.
#'@param end end date of records for magnetic activities.See examples
#'@param n.stations number of stations included in the study
#'@param xlab x-axis label
#'@param ylab y-axis label
#'@param ... additional arguments
#'
#'@details This is an interface for visualizing WISA with ggplot2.
#'\code{SIplot} uses \code{as.POSIXct} to convert the index object into a data
#'frame with datetime. It essentially uses the mapping \code{geom_line} and adds
#'aes(x, y) to visualize the series.
#'
#'@seealso \code{\link{SAIndex}}, \code{\link{preindexplot}}
#'
#'@examples
#'\dontrun{
#'## sample records for one week period
#'data(record)
#'coord=matrix(c(124.43, 19.23, 53.77,140.18,68.68,202.00,71.89,293.85),nrow=2,ncol=4)
#'
#'index.sample<- SAIndex(record, coord, wf="la8")
#'
#'## example dates
#'start.date="2001-3-1"
#'end.date="2001-4-30"
#'n.station=4
#'station.names=c("HER","KAK","HON","SJG")
#'
#'## plot SI
#'SIplot(index.sample$SI, type=1, start=start.date, end=end.date, 
#'main="WISA estimation", xlab="Datetime", ylab="iWISA Estimation")
#'}
#'@import ggplot2
#'@export
SIplot <- function(x, type=NULL, main = NULL, xlab = NULL,
                         ylab = NULL, start = NULL,end=NULL,n.stations=NULL,...)
{
#index.si<-data[1]
#build the data frame

days<-seq(as.Date(start), as.Date(end), by = "day")
date.frame.si<- data.frame(Date=as.Date(rep(days,each=1440)),
               Time=format(as.difftime(0:1439, units = "mins") +
                                 as.POSIXct(start), "%H:%M"))
index.value.si<-cbind(date.frame.si,x)
index.value.si$datetime <- as.POSIXct(paste(index.value.si$Date,index.value.si$Time),tz="UCT")

si<-index.value.si[3:4]
#list.files(system.file('enc', package = 'grDevices'))

gg <- ggplot(si)
gg <- gg + geom_line(lwd = type,aes_string(x='datetime', y='x', color=24))
gg <- gg + labs(x=xlab, y=ylab, title=main)
gg <- gg + theme_bw(base_family="Helvetica")
gg <- gg + theme(axis.ticks.y=element_blank())
gg <- gg + theme(panel.border=element_blank())
gg <- gg + theme(legend.key=element_blank())
gg <- gg + theme(legend.position="none")
gg

}
