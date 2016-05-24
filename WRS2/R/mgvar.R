mgvar<-function(m,se=FALSE,op=0,cov.fun=covmve,SEED=TRUE){
#
# Find the center of a scatterplot, add point that
# increases the generalized variance by smallest amount
# continue for all points
# return the generalized variance
#  values corresponding to each point.
# The central values and point(s) closest to it get NA
#
# op=0 find central points using pairwise differences
# op!=0 find central points using measure of location
# used by cov.fun
#
# choices for cov.fun include
# covmve
# covmcd
# tbs (Rocke's measures of location
# rmba (Olive's median ball algorithm)
#
if(op==0)temp<-apgdis(m,se=se)$distance
if(op!=0)temp<-out(m,cov.fun=cov.fun,plotit=FALSE,SEED=SEED)$dis
flag<-(temp!=min(temp))
temp2<-temp
temp2[!flag]<-max(temp)
flag2<-(temp2!=min(temp2))
flag[!flag2]<-F
varvec<-NA
while(sum(flag)>0){
ic<-0
chk<-NA
remi<-NA
for(i in 1:nrow(m)){
if(flag[i]){
ic<-ic+1
chk[ic]<-gvar(rbind(m[!flag,],m[i,]))
remi[ic]<-i
}}
sor<-order(chk)
k<-remi[sor[1]]
varvec[k]<-chk[sor[1]]
flag[k]<-F
}
varvec
}
