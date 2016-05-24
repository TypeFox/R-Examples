get.peaks<-function(data,channel,trials=NULL,
trial.cn="Trial",time.cn="Time",add.lines=NULL,...)
{
	if(is.null(trials)){
		trials<-sort(unique(data[,trial.cn]))
	}		
	xy<-list()
	for(tr in trials){
		cat("fetching peaks for trial",tr,"...\n")
		tmp<-data[data[,trial.cn]==tr,]
		plot(tmp[,time.cn],tmp[,channel],type="l",main=paste("Trial",tr),
			xlab="Time",ylab="Amplitude",...)
		if(!is.null(add.lines)){
			for(al in 1:length(add.lines)){
				lines(tmp[,time.cn],tmp[,add.lines[[al]][1]],
					col=add.lines[[al]][2],...)
			}
		}
		xy.tmp<-locator(...)
		if(!is.null(xy.tmp)){
		time<-tmp[,time.cn]
		MYTIME<-vector("numeric")
		for(k in 1:length(xy.tmp$x)){
			mytime<-time[time<=xy.tmp$x[k]]
			mytime<-mytime[length(mytime)]
			mytime<-c(mytime,time[time>xy.tmp$x[k]][1])
			mydist<-abs(abs(xy.tmp$x[k])-abs(mytime[1]))
			mydist<-c(mydist,abs(abs(xy.tmp$x[k])-abs(mytime[2])))
			MYTIME<-c(MYTIME,mytime[which(mydist==min(mydist))])
		}
		}else{
			cat("    no xy-coordinates fetched for trial",tr,"\n")
			MYTIME<-NA
		}
		xy[[which(trials==tr)]]<-MYTIME
	}
	names(xy)<-paste("trial",trials,sep="")
	dev.off()
	return(invisible(xy))
}

