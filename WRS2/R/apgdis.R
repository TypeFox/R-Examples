apgdis<-function(m,est=sum,se=TRUE,...){
#
# For multivariate data,
# compute distance between each pair
# of points and measure depth of a point
# in terms of its  distance to all
# other points
#
#  Using se=T ensures that ordering of distance
# will not change with a change in scale.
#
#  m is an n by p matrix
#
m<-elimna(m)  # eliminate any missing values
temp<-0
if(se){
for(j in 1:ncol(m))m[,j]<-(m[,j]-median(m[,j]))/mad(m[,j])
}
for(j in 1:ncol(m)){
disx<-outer(m[,j],m[,j],"-")
temp<-temp+disx^2
}
temp<-sqrt(temp)
dis<-apply(temp,1,est,...)
temp2<-order(dis)
center<-m[temp2[1],]
list(center=center,distance=dis)
}

