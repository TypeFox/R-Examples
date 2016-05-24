dmean<-function(m,tr=.2,dop=1,cop=2){
#
# Compute multivariate measure of location
# using Donoho-Gasko method.
#
# dop=1, use fdepth to compute depths
# dop=2, use fdepthv2  to compute depths
#
# cop=1, Tukey median; can't be used here.
# cop=2, use MCD in fdepth
# cop=3, use marginal medians in fdepth
# cop=4, use MVE in fdepth
#
if(is.list(m))m<-matl(m)
if(!is.matrix(m))stop("Data must be stored in a matrix or in list mode.")
if(ncol(m)==1){
if(tr==.5)val<-median(m)
if(tr>.5)stop("Amount of trimming must be at most .5")
if(tr<.5)val<-mean(m,tr)
}
if(ncol(m)>1){
temp<-NA
if(ncol(m)!=2){
# Use approximate depth
if(dop==1)temp<-fdepth(m,plotit=FALSE,cop=cop)
if(dop==2)temp<-fdepthv2(m)
}
#  Use exact depth if ncol=2
if(ncol(m)==2){
for(i in 1:nrow(m))
temp[i]<-depth(m[i,1],m[i,2],m)
}
mdep<-max(temp)
flag<-(temp==mdep)
if(tr==.5){
if(sum(flag)==1)val<-m[flag,]
if(sum(flag)>1)val<-apply(m[flag,],2,mean)
}
if(tr<.5){
flag2<-(temp>=tr)
if(sum(flag2)==0 && sum(flag)>1)val<-apply(as.matrix(m[flag,]),2,mean)
if(sum(flag2)==0 && sum(flag)==1)val=m[flag,]
if(sum(flag2)==1)val<-m[flag2,]
if(sum(flag2)>1)val<-apply(m[flag2,],2,mean)
}}
val
}
