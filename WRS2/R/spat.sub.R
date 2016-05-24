spat.sub<-function(x,theta){
xx<-x
for(i in 1:ncol(x))xx[,i]<-x[,i]-theta[i]
xx<-xx^2
temp<-sqrt(apply(xx,1,sum))
val<-mean(temp)
val
}