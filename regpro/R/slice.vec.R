slice.vec<-function(x,coordi=1,p=0.5)
{
d<-dim(x)[2]
center<-matrix(0,d,1)
for (i in 1:d) center[i]<-quantile(x[,i],p=p)  #center<-colMeans(x)
direc<-rep(0,d)
direc[coordi]<-1
radius<-(max(x[,coordi])-min(x[,coordi]))/2

return(list(center=center,direc=direc,radius=radius))
}



