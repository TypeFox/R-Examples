#'Estimate Index and Preindex of Storm Activity
#'
#'Function to estimate index of global storm activity, WISA (wavelet index of storm activity), and values of storm activity preindex of individual stations Hermanus (HER), station KAkioka (KAK), station Honolulu (HON) and station San Juan (SJG). The WISA can be constructed practically over any period of time, but to obtain values similar to the standard Dst, at least 2 months worth of data are needed.
#'
#'@param data a matrix or a data frame of H-component values. Each column should contain records from an individual station.
#'@param coord a matrix of colatitude and longitude of each station used in a study
#'@param wf type of wavelet filter used to separate magnetic storm activity from the solar quiet variation. LA(8) filter is the default.
#'
#'@details \code{SAIndex} is a function for estimating the index of magnetic storm activity. The function returns the index of global storm activity as well as the preindex for each of the individual stations used.
#'@return SI the index of global storm activity, WISA
#'@return p.SQ the preindex of the storm activity constructed using wavelet-based MRA. It is the contribution of a station to the storm signature.
#'@references Jach,A., P.Kokoszka, J.Sojka, and L.Zhu(2006),Wavelet-based index of magnetic storm activity, J. Geophys. Res.,111, A09215, doi:10.1029/2006JA011635.
#'@references Maslova, I.,P.Kokoszka, J.Sojka, and L. Zhu(2009),Removal of nonconstant daily variation by means of wavelet and functional data analysis, J. Geophys. Res.,114, A03202, doi:10.1029/2008JA013685.
#'@author Inga Maslova
#'@author Wen Xiao
#'@seealso \code{\link{magnetic.latitude}}, \code{\link{SIplot}}, \code{\link{preindexplot}}
#'
#'@examples
#'\dontrun{
#'## data sets.
#'data(record)
#'## Colatitude and longitude of geomagnetic observatories.
#'coord=matrix(c(124.43, 19.23, 53.77,140.18,68.68,202.00,71.89,293.85),nrow=2,ncol=4)
#'
#'
#'## estimate WISA and preindex for each stations.
#'index.sample<- SAIndex(record, coord, wf="la8")
#'
#'## plot SI
#'start.date="2001-3-1"
#'end.date="2001-4-30"
#'n.station=4
#'station.names=c("HER","KAK","HON","SJG")
#'
#'
#'SIplot(index.sample$SI, type=1,start=start.date, end=end.date, 
#'main="WISA estimation", xlab="Datetime", ylab="iWISA Estimation")
#'}
#'
#'@export


SAIndex <-
function(data, coord, wf="la8"){
###########################################
# filter each Dst station and adjust for the magnetic latitude and
# subtract the average
###########################################
J0<-7
J1<-10
boundary<-"reflection"
periods<-60*c(8,12,24)
quantile<-0.98
n.levels<-11

SAI<-WISA(data, wf1 = "la8", n.levels = n.levels, J0 = J0, J1 = J1,
            boundary = boundary, quantile = quantile)

n.station<-dim(data)[2] #each column of the data matrix are records for each individual station
N<-dim(data)[1]         # number of observations for each station

index<-matrix(0, ncol=n.station, nrow=N)
final.index<-matrix(0, ncol=n.station, nrow=N)

# make latitude adjustment.
for(i in 1:n.station){
# adjust for the magnetic latitude
index[,i]<-SAI$preindex[,i]/cos(magnetic.latitude(coord[1,i],coord[2,i])$Phi*2*pi/360)
#}

# center: subtract the mean
for(i in 1:n.station){
# adjust for the magnetic latitude
final.index[,i]<-index[,i]-mean(index[,i])
#}
###########################################

SI<-numeric(N)
SI<-apply(final.index, 1, mean)
}
}
list(SI=SI, p.SQ=SAI$pseudo.sq)
#multi_return <- function() {
#   my_list <- list(SI=SI, p.SQ=SAI$pseudo.sq)
#   return(my_list)
#}
#multi_return()
}
