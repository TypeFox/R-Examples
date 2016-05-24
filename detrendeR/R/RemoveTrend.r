RemoveTrend = function (rw, method=c("Spline", "Spline%", "ModNegExp", "Mean"),BandwidthPerc=2/3, Bandwidth=32, P=1/2){

BandwidthPerc = as.numeric(BandwidthPerc)

if(method=="Spline") {curve<-SPLINE(rw, bandwidth=Bandwidth, p=P); return(curve)}
if(method=="Spline%") {curve<-SPLINE(rw, bandwidth=length(rw)*BandwidthPerc, p=P); return(curve)}
if(method=="ModNegExp") {curve<-NegExp(rw); return(curve)}
if(method=="Mean") {series<-rw[!is.na(rw)]; curve<-rep(mean(series),length(series)) ; rw[!is.na(rw)]=curve; return(rw)}

}
