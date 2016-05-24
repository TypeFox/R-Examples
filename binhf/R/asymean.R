"asymean" <-
function(xgrid=seq(0,1,length=21),ygrid=seq(0,1,length=21),binsize=32){

zetam1m2<-matrix(0,length(xgrid),length(ygrid))

for (i in 1:length(xgrid)){
for (j in 1:length(xgrid)){
zetam1m2[i,j]<-(ygrid[j]-xgrid[i])/sqrt((ygrid[j]+xgrid[i])*(2-(ygrid[j]+xgrid[i]))/(2*binsize))
}
}
zetam1m2[which(abs(zetam1m2)==Inf)]<-0
zetam1m2[which(is.na(zetam1m2))]<-0

zetam1m2
}

