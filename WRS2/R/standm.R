standm<-function(x,locfun=lloc,est=mean,scat=var,...){
# standardize a matrix x
#
x=elimna(x)
x=as.matrix(x)
m1=lloc(x,est=est)
v1=apply(x,2,scat)
p=ncol(x)
for(j in 1:p)x[,j]=(x[,j]-m1[j])/sqrt(v1[j])
x
}
