mom<-function(x,bend=2.24,na.rm=TRUE){
#
#  Compute MOM-estimator of location.
#  The default bending constant is 2.24
#
if(na.rm)x<-x[!is.na(x)] #Remove missing values
flag1<-(x>median(x)+bend*mad(x))
flag2<-(x<median(x)-bend*mad(x))
flag<-rep(T,length(x))
flag[flag1]<-F
flag[flag2]<-F
mom<-mean(x[flag])
mom
}