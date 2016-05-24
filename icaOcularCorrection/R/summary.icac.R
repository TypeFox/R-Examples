summary.icac<-function(object,noise.sig=NULL,ic=NULL,print=TRUE,...){
	if(!"icac"%in%class(object))stop("object not of class \"icac\"\n")
	if(missing(object))stop("please supply an \"icac\" object\n")
	#if(missing(noise.sig))stop("please supply a noise signal\n")
	ci<-object$correction.info
	if(!is.null(noise.sig)){
		ci<-ci[ci$NoiseSignal==noise.sig,]
	}
	if(!is.null(ic)){
		ci<-ci[ci$IC==ic,]
		if(is.null(noise.sig)){
			ns<-sort(unique(ci$NoiseSignal))
			for(kk in ns){
				if(which(ns==kk)==1){
					smry<-data.frame(IC=ic,NoiseSignal=kk,
						NumTrials=nrow(ci[ci$NoiseSignal==kk,]),
						MeanCorr=mean(abs(ci[ci$NoiseSignal==kk,]$Corr),
						na.rm=TRUE))
				}else{
					smry<-rbind(smry,data.frame(IC=ic,NoiseSignal=kk,
						NumTrials=nrow(ci[ci$NoiseSignal==kk,]),
						MeanCorr=mean(abs(ci[ci$NoiseSignal==kk,]$Corr),
						na.rm=TRUE)))
				}
			}
			if(print){
			mytitle<-"SUMMARY FOR"
				mytitle<-paste(mytitle,"IC =",ic)
				
				cat("\n")
				cat(paste(mytitle,":\n",sep=""))
				print(smry)
				cat("----------\n")
				cat("NOTE: More infor available in ICAC_OBJECT$correction.info\n")
			}
			return(invisible(smry))
		}
	} 

	tmp<-ci

	ics.pres<-sort(unique(tmp$IC))
	for(ii in ics.pres){
		if(which(ics.pres==ii)==1){
			smry<-data.frame(IC=ii,NumTrials=nrow(tmp[tmp$IC==ii,]),
				MeanCorr=mean(abs(tmp[tmp$IC==ii,]$Corr),na.rm=TRUE))
		}else{
			smry<-rbind(smry,data.frame(IC=ii,
				NumTrials=nrow(tmp[tmp$IC==ii,]),
				MeanCorr=mean(abs(tmp[tmp$IC==ii,]$Corr),na.rm=TRUE)))
		}
	}
	smry<-smry[order(smry$NumTrials,smry$MeanCorr,decreasing=TRUE),]
	rownames(smry)<-1:nrow(smry)

	if(print){
		mytitle<-"SUMMARY"
		if(!is.null(noise.sig)){
			mytitle<-paste(mytitle,"FOR NOISE SIGNAL =",noise.sig)
		}
		if(!is.null(ic)){
			mytitle<-paste(mytitle,"AND IC =",ic)
		}
		
		cat("\n")
		cat(paste(mytitle,":\n",sep=""))
		print(smry)
		cat("----------\n")
		cat("NOTE: More infor available in ICAC_OBJECT$correction.info\n")
	}
	
	return(invisible(smry))
}
