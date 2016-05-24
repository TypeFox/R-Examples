pdis<-function(m,MM=FALSE,cop=3,dop=1,center=NA){
#
# Compute projection distances for points in m
#
#
#
#  MM=F  Projected distance scaled
#  using interquatile range.
#  MM=T  Scale projected distances using MAD.
#
#  There are five options for computing the center of the
#  cloud of points when computing projections:
#  cop=1 uses Donoho-Gasko median
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#  cop=5 uses skipped mean
#
m<-elimna(m) # Remove missing values
m<-as.matrix(m)
if(ncol(m)==1){
if(is.na(center[1]))center<-median(m)
dis<-abs(m[,1]-center)
if(!MM){
temp<-idealf(dis)
pdis<-dis/(temp$qu-temp$ql)
}
if(MM)pdis<-dis/mad(dis)
}
if(ncol(m)>1){
if(is.na(center[1])){
if(cop==1)center<-dmean(m,tr=.5,dop=dop)
if(cop==2)center<-cov.mcd(m)$center
if(cop==3)center<-apply(m,2,median)
if(cop==4)center<-cov.mve(m)$center
if(cop==5)center<-smean(m)
}
dmat<-matrix(NA,ncol=nrow(m),nrow=nrow(m))
for (i in 1:nrow(m)){
B<-m[i,]-center
dis<-NA
BB<-B^2
bot<-sum(BB)
if(bot!=0){
for (j in 1:nrow(m)){
A<-m[j,]-center
temp<-sum(A*B)*B/bot
dis[j]<-sqrt(sum(temp^2))
}
if(!MM){
temp<-idealf(dis)
dmat[,i]<-dis/(temp$qu-temp$ql)
}
if(MM)dmat[,i]<-dis/mad(dis)
}}
pdis<-apply(dmat,1,max,na.rm=TRUE)
}
pdis
}
