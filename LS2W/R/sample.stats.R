`sample.stats`<-function (x, n=25, size=64)
{
dim1<-nrow(x)
dim2<-ncol(x)
xcoord<-sample.int(dim1-size, n)
ycoord <- sample.int(dim2-size,n)
J<-log(size,2)
feature.matrix<-matrix(0,nrow=n, ncol=3*J)
for (i in 1:n){
	sam <- x[xcoord[i]:(xcoord[i]+size-1), ycoord[i]:(ycoord[i]+size-1)]
        sam.ls2w<-cddews(sam, filter.number=1, family="DaubExPhase", levels=3:5)
        feature.matrix[i,]<-apply(sam.ls2w$S,1,sum)
}
return(feature.matrix)
}
