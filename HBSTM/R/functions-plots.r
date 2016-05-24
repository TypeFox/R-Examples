

##################################################################################################################
######################################      Plot functions     ###################################################
##################################################################################################################

.HBSTM.resid.plot=function(object,point,ARlags,ARperiod){
	if(missing(point)) point=sample(1:dim(object@newGrid)[1],1)	
	if(is.na(object@fitted[1,1,1]))	stop("The object does not have stored fitted values\n")
	Et=object@Zt-object@K%*%object@fitted[,,1]+object@Parameters@sigma2E
	Et=Et[point,]
	dev.new()
	layout.show(layout(mat=matrix(c(1,2,3,3),ncol=2,nrow=2,byrow=TRUE)))
	hist(Et,main=paste("Residuals histogram in s=",point,sep=""))
	plot(Et,main=paste("Residuals in s=",point,sep=""))
	abline(h=0,col=2)
	qqnorm(Et,main=paste("Normal Q-Q Plot Residuals in s=",point,sep=""))
	qqline(Et,col=2,lwd=2)	
	dev.new()
	if(missing(ARlags)) ARlags=NULL
	if(missing(ARperiod)){
		color=1
	}else{
		color=c(rep(1,ARperiod-1),2)
	}
	par(mfrow=c(2,2))
	acf(Et,lag.max=ARlags,col=color,main=paste("ACF of the residuals in s=",point,sep=""))
	pacf(Et,lag.max=ARlags,col=color,main=paste("ACF of the residuals in s=",point,sep=""))
	acf(Et^2,lag.max=ARlags,col=color,main=paste("ACF of the residuals^2 in s=",point,sep=""))
	pacf(Et^2,lag.max=ARlags,col=color,main=paste("ACF of the residuals^2 in s=",point,sep=""))	

	return(invisible())
}
setMethod(f="plotRes",signature=c("HBSTM","ANY","ANY","ANY"),definition=.HBSTM.resid.plot)



plotCI.aux=function(name,spatTemp,mat,labnames,las){
	if(missing(las)) las=0
	plot(x=1:length(spatTemp),ylim=c(min(mat[,2]),max(mat[,3])),xlim=c(0.5,length(spatTemp)+0.5),type="n",ylab="Values",xlab="",main=name,xaxt="n")
	axis(1,at=1:length(spatTemp),labels=labnames,las=las)
	for(i in 1:length(spatTemp)){
		arrows(x0=i,y0=mat[i,2],x1=i,y1=mat[i,3],code=3,angle=90,length=0.1)
		points(x=i,y=mat[i,1])
	}
}	



plotCI.spatial=function(name,spat,labs,inter,report=FALSE,nplot){
	if(length(spat)<5){
		if(report){
			jpeg(paste("plot-",nplot,".jpeg",sep=""))
		}else{
			dev.new()
		}
		par(omi=c(0,0,0.2,0))
		if(length(spat)==1)	layout(matrix(1,nrow=1,ncol=1))
		if(length(spat)==2)	layout(matrix(1:2,nrow=2,ncol=1))
		if(length(spat)==3)	layout(matrix(c(1:3,0),nrow=2,ncol=2))
		if(length(spat)==4)	layout(matrix(c(1:4),nrow=2,ncol=2))
	}else{
		if(length(spat)<7){
			if(report){
				jpeg(paste("plot-",nplot,".jpeg",sep=""))
			}else{
				dev.new(width=15)
			}
			par(omi=c(0,0,0.2,0))
			layout(matrix(c(1:length(spat)),nrow=2,ncol=3))
		}else{
			if(report){
				jpeg(paste("plot-",nplot,".jpeg",sep=""))
			}else{
				dev.new(width=20)
			}
			layout(matrix(c(1:length(spat)),nrow=2,ncol=4))
		}
	}	
	for(k in 1:length(spat)){
		mat=as.data.frame(matrix(NA,ncol=3,nrow=dim(spat[[k]])[2]))
		for(i in 1:dim(spat[[k]])[2]){
			mat[i,]=mat.aux(spat[[k]][,i],inter)
		}
		plotCI.aux(name=names(spat)[k],spatTemp=1:dim(spat[[k]])[2],mat=mat,labnames=paste(names(spat)[k],1:dim(spat[[k]])[2],sep=""))
	}
	if(length(spat)<5){
		mtext(name,line=-1,outer=TRUE,font=2,cex=1.2)
	}else{
		mtext(name,line=-1.5,outer=TRUE,font=2,cex=1.2)
	}
}

auxplotFit=function(long,lat,Zt,main){
	surf=linearInterp(long ,lat,Zt,gridPoints=as.integer(length(Zt)*0.8))
	image(surf$x,surf$y,surf$z, add = F,col=topo.colors(12),xlab="LON",ylab="LAT",main=main)
	contour(surf$x,surf$y,surf$z,add = TRUE)		
}

.HBSTM.plotFit=function(object,time,values){
	
	if(missing(time)) time=dim(object@Zt)[2]
	if(missing(values)) values=FALSE
	
	longlat=object@newGrid
	Zt=object@Zt[,time]

	dev.new()
	layout(matrix(c(1,1,2,2,0,3,3,0),nrow=2,ncol=4,byrow=TRUE))
	auxplotFit(long=longlat[,1],lat=longlat[,2],Zt=object@Zt[,time],main=paste("Zt observations in t =",time))
	
	auxplotFit(long=longlat[,1],lat=longlat[,2],Zt=object@K%*%object@fitted[,time,1]+object@Parameters@sigma2E,main=paste("Zt estimations in t =",time))
		
	auxplotFit(long=longlat[,1],lat=longlat[,2],Zt=object@Zt[,time]-(object@K%*%object@fitted[,time,1]+object@Parameters@sigma2E),main=paste("Residuals in t =",time))
	
	if(values){
		return(data.frame(Zt=object@Zt[,time],EstZt=object@K%*%object@fitted[,time,1]+object@Parameters@sigma2E,Et=object@Zt[,time]-object@K%*%object@fitted[,time,1]+object@Parameters@sigma2E))
	}else{
		return(invisible())
	}
	
}
setMethod(f="plotFit",signature=c("HBSTM","ANY","ANY"),definition=.HBSTM.plotFit)