readGDALbbox<-function(gdal,spo,mar=2,...){
	mybb<-bbox(spo)
	x<-GDALinfo(gdal)
	
	myx<-seq(x[[4]],x[[4]]+x[[2]]*x[[6]],by=x[[6]])
	myy<-seq(x[[5]]+x[[1]]*x[[7]],x[[5]],by=-x[[7]])

	xcol<-which(myx>=mybb[1,1] & myx<=mybb[1,2])
	ycol<-which(myy>=mybb[2,1] & myy<=mybb[2,2])

  x<-readGDAL(gdal,offset=c(min(ycol)-mar,min(xcol)-mar),region.dim=c(length(ycol)+mar,length(xcol)+mar),...)
  return(x)
}

