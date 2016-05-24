fdepthv2<-function(m,pts=NA,plotit=TRUE){
#
# Determine depth of points in pts relative to
# points in m
#
# Draw a line between each pair of distinct points
# and determine depth of the projected points.
# The final depth of a point is its minimum depth
# among all projections.
#
# This function is slower than fdepth and requires
# space for a nc by nc matrix, nc=(n^2-n)/2.
# But it allows
# data to have a singular covariance matrix
# and it provides a more accurate approximation of
# halfspace depth.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data and pts=NA
#
#  When plotting,
#  center is marked with a cross, +.
#
m<-elimna(m) # Remove missing values
if(!is.na(pts[1]))remm<-m
if(!is.matrix(m))dep<-unidepth(m)
if(is.matrix(m)){
nm<-nrow(m)
nt<-nm
nm1<-nm+1
if(!is.na(pts[1])){
if(ncol(m)!=ncol(pts))stop("Number of columns of m is not equal to number of columns for pts")
nt<-nm+nrow(pts)
}}
if(ncol(m)==1)depth<-unidepth(m)
if(ncol(m)>1){
m<-elimna(m) # Remove missing values
nc<-(nrow(m)^2-nrow(m))/2
if(is.na(pts[1]))mdep <- matrix(0,nrow=nc,ncol=nrow(m))
if(!is.na(pts[1])){
mdep <- matrix(0,nrow=nc,ncol=nrow(pts))
}
ic<-0
for (iall in 1:nm){
for (i in 1:nm){
if(iall < i){
ic<-ic+1
B<-m[i,]-m[iall,]
dis<-NA
BB<-B^2
bot<-sum(BB)
if(bot!=0){
if(is.na(pts[1])){
for (j in 1:nrow(m)){
A<-m[j,]-m[iall,]
temp<-sum(A*B)*B/bot
dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
}}
if(!is.na(pts[1])){
m<-rbind(remm,pts)
for (j in 1:nrow(m)){
A<-m[j,]-m[iall,]
temp<-sum(A*B)*B/bot
dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
}}
#
# For ic_th projection, store depths of
# points in mdep[ic,]
#
if(is.na(pts[1]))mdep[ic,]<-unidepth(dis)
if(!is.na(pts[1])){
mdep[ic,]<-unidepth(dis[1:nm],dis[nm1:nrow(m)])
}}
if(bot==0)mdep[ic,]<-rep(0,ncol(mdep))
}}}
dep<-apply(mdep,2,min)
}
if(ncol(m)==2 &&is.na(pts[1])){
flag<-chull(m)
dep[flag]<-min(dep)
}
if(ncol(m)==2){
if(is.na(pts[1]) && plotit){
plot(m)
x<-m
temp<-dep
flag<-(temp>=median(temp))
xx<-x[flag,]
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
xord<-order(xx[,1])
xx<-xx[xord,]
temp<-chull(xx)
lines(xx[temp,])
lines(xx[c(temp[1],temp[length(temp)]),])
}}
dep
}
