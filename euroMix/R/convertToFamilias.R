convertToFamilias=function(infile,outfile=paste("out",infile,sep="")){
dat1=read.table(file=infile,header=T)
markerNames=colnames(dat1)[-(1:2)]
dimDat=dim(dat1)
d1=dimDat[1]
d2=dimDat[2]
seqOdd=seq(1,d1,by=2)
seqEven=seq(2,d1,by=2)
dat2=cbind(dat1[seqOdd,],dat1[seqEven,])
d3=dim(dat2)[2]
index=1:d3
index[seq(1,d3,by=2)]=1:d2
index[seq(2,d3,by=2)]=(1:d2)+d2
dat2=dat2[,index]
dat2=dat2[,-2]
markerNames=rep(markerNames,each=2)
markerNames=paste(markerNames,1:2)
line1=c("Name","Amel1","Amel2",markerNames)
sink(outfile)
cat(line1,"\n",sep="\t")
write.table(dat2,file=outfile,row.names=FALSE,quote=FALSE,sep="\t",append=TRUE,col.names=FALSE)
sink()
cat(outfile," printed","\n")
}


