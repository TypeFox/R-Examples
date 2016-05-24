lagSelect<-function(sampleData,sampleLocation=NULL,sampleTime=NULL,LocationNames=NULL,lagmin,lagmax,by,
	option=1,numTests=99,plot=T,anisotropic=F,...){
	len<-seq(lagmin,lagmax,by)
	if(option %in% c(2,4)) maxdist<-max(vegdist(sampleTime,method='eu'))	
	   else maxdist<-max(vegdist(sampleLocation,method='eu'))
	tests<-vector("list",length=length(len))
	for (i in c(1:length(len))){
		tests[[i]]<-list(lagSize=len[i],
			community.correlogram=commcorrelogram(sampleData=sampleData
			,sampleLocation=sampleLocation,sampleTime=sampleTime
			,LocationNames=LocationNames,option=option
			,lagNumber=ceiling(maxdist/(2*len[i]))
			,lagSize=len[i],lagTol=len[i]/2
			,numTests=numTests,anisotropic=anisotropic,...))
	}
	if(plot==TRUE)
	for (i in c(1:length(len))){
		dev.new()
		par(mfrow=c(2,1),mar=c(4,4,3,2))
		plot(tests[[i]][[2]]@community.correlogram[,c(1,3)],xlab='lag distance'
		,ylab=colnames(tests[[i]][[2]]@community.correlogram)[3]
		,main=paste('Lag Size',tests[[i]][[1]],sep='='),type='o'
		,cex.axis=0.7,pch=as.numeric(tests[[i]][[2]]@community.correlogram[,2]<30)*7+1)
		par(mar=c(4,4,3,2))
		plot(tests[[i]][[2]]@community.correlogram[,c(1,4)],xlab='lag distance'
		,ylab='significance',main='',type='o',cex.axis=0.7
		,pch=as.numeric(tests[[i]][[2]]@community.correlogram[,2]<30)*7+1)
	}
return(tests)
}
