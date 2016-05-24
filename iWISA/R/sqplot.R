#'Plot of Non-constant Solar Quiet Daily Variation
#'
#'This function visualizes the estimated non-constant solar quiet daily variation of each station.
#'Users would define how many graphs to be put into one page, The default set is 2 graphs per page.
#'
#'@param x estimation of solar quite variations from \code{SQ}
#'@param start start date of records for magnetic activities. See examples
#'@param end end date of records for magnetic activities. See examples
#'@param n.station number of studied stations
#'@param graphs.per.page how many graphs combined in one plot page. Default number is 2
#'@param station.names NULL (indicating no station names) or a vector strings for station names.
#'@param Title title of sqplot
#'@param ... additional arguments
#'
#'@details This function is used to visualize estimated Sq variation. Specifying the names of each stations would
#'add station names at y-axis.
#'
#'@seealso \code{\link{SQ}}
#'@examples
#'## data sets.
#'
#'## Colatitude and longitude of geomagnetic observatories.
#'\dontrun{
#'coord=matrix(c(124.43, 19.23, 53.77,140.18,68.68,202.00,71.89,293.85),nrow=2,ncol=4)
#'
#'
#'## Estimation of SI index
#'index.sample<- SAIndex(record, coord, wf="la8")
#'si.v<-index.sample$SI
#'
#'
#'## generate datetime for one week period
#'start.date="2001-3-1"
#'end.date="2001-4-30"
#'
#'## estimation of sq
#'Sq<- SQ (datasq, si.v=si.v, wf = "la8")
#'
#'sqplot(Sq, Title="Sq variation", start = start.date, end=end.date, n.station=4,
#'       graphs.per.page=2, station.names=c("ABG","PHU","TUC","FRD"))
#'}
#'@export
#'@importFrom graphics axis.POSIXct par plot 
#'
sqplot <- function(x, start = NULL,end=NULL,n.station=NULL,
graphs.per.page=2, station.names=NULL, Title=NULL,...)
{
#data(datasq, package="iWISA")
#N<-nrow(data.frame(datasq))
N<-nrow(data.frame(x))
index.sq<-matrix(0, ncol=n.station, nrow=N)
colnames(index.sq) <-station.names
# paste("Station", station.names, sep = " ")
for (i in 1:n.station){
#index.sq[,i]<-datasq[,i]}
index.sq[,i]<-x[,i]}
days<-seq(as.Date(start), as.Date(end), by = "day")
dates.value<- data.frame(Date=as.Date(rep(days,each=1440)),
               Time=format(as.difftime(0:1439, units = "mins") +
                                 as.POSIXct(start), "%H:%M"))
index.value<-cbind(dates.value, index.sq)
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
