onestep<-function(x,bend=1.28,na.rm=FALSE,MED=TRUE){
#
#  Compute one-step M-estimator of location using Huber's Psi.
#  The default bending constant is 1.28
#  
#  MED=TRUE: initial estimate is the median
#  Otherwise use modified one-step M-estimator
#
if(na.rm)x<-x[!is.na(x)]
if(MED)init.loc=median(x)
if(!MED)init.loc=mom(x,bend=bend)
y<-(x-init.loc)/mad(x)  #mad in splus is madn in the book.
A<-sum(hpsi(y,bend))
B<-length(x[abs(y)<=bend])
onestep<-median(x)+mad(x)*A/B
onestep
}
