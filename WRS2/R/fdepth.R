fdepth<-function(m,pts=NA,plotit=TRUE,cop=3,center=NA,xlab="VAR 1",
ylab="VAR 2"){
#
# Determine depth of points in pts,  relative to
# points in m. If pts is not specified,
# depth of all points in m are determined.
#
# m and pts can be vectors or matrices with
# p columns (the number of variables).
#
# Determine center, for each point, draw a line
# connecting it with center, project points onto this line
# and determine depth of the projected points.
# The final depth of a point is its minimum depth
# among all projections.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data and pts=NA
#
#  There are three options for computing the center of the
#  cloud of points when computing projections, assuming center=NA:
#
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#
#  If a value for center is passed to this function,
#  this value is used to determine depths.
#
#  When plotting,
#  center is marked with a cross, +.
#

if(cop!=2 && cop!=3 && cop!=4)stop("Only cop=2, 3 or 4 is allowed")
if(is.list(m))stop("Store data in a matrix; might use function listm")
m<-as.matrix(m)
pts<-as.matrix(pts)
if(!is.na(pts[1]))remm<-m
nm<-nrow(m)
nm1<-nm+1
if(!is.na(pts[1])){
if(ncol(m)!=ncol(pts))stop("Number of columns of m is not equal to number of columns for pts")
}
m<-elimna(m) # Remove missing values
m<-as.matrix(m)
if(ncol(m)==1)dep<-unidepth(as.vector(m[,1]),pts=pts)
if(ncol(m)>1){
if(is.na(center[1])){
if(cop==2){
center<-cov.mcd(m)$center
}
if(cop==4){
center<-cov.mve(m)$center
}
if(cop==3){
center<-apply(m,2,median)
}}
if(is.na(pts[1])){
mdep <- matrix(NA,nrow=nrow(m),ncol=nrow(m))
}
if(!is.na(pts[1])){
mdep <- matrix(NA,nrow=nrow(m),ncol=nrow(pts))
}
for (i in 1:nrow(m)){
B<-m[i,]-center
dis<-NA
BB<-B^2
bot<-sum(BB)
if(bot!=0){
if(is.na(pts[1])){
for (j in 1:nrow(m)){
A<-m[j,]-center
temp<-sum(A*B)*B/bot
dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
}}
if(!is.na(pts[1])){
m<-rbind(remm,pts)
for (j in 1:nrow(m)){
A<-m[j,]-center
temp<-sum(A*B)*B/bot
dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
}}
#
# For ith projection, store depths of
# points in mdep[i,]
#
if(is.na(pts[1]))mdep[i,]<-unidepth(dis)
if(!is.na(pts[1])){
mdep[i,]<-unidepth(dis[1:nm],dis[nm1:nrow(m)])
}}
if(bot==0)mdep[i,]<-rep(0,ncol(mdep))
}
dep<-apply(mdep,2,min)
if(ncol(m)==2 && is.na(pts[1])){
flag<-chull(m)
dep[flag]<-min(dep)
}
}
if(ncol(m)==2){
if(is.na(pts[1]) && plotit){
plot(m,xlab=xlab,ylab=ylab)
points(center[1],center[2],pch="+")
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
dep<-round(dep*nrow(m))/nrow(m)
dep
}
