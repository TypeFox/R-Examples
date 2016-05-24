plot_nic<-function(x,
		   data,
		   ic,
		   trial,
		   noise.sig,
		   threshold=x$threshold,
		   col=c("black","blue"),
		   main=NULL,
		   xlab=NULL,
		   ylab=NULL,
           xlim=NULL,
		   ylim=NULL,
		   cex=0.7,
		   trial.cn="Trial",
		   ...
){
	if(missing(x))stop("Please provide a object \"x\" of class \"icac\".\n")
	if(missing(data))stop("Please provide an object \"data\" containing the original uncorrected data.\n")
	if(missing(ic))stop("Please provide an IC.\n")
	if(missing(trial))stop("Please provide trial.\n")
	if(missing(noise.sig)){
    	if(length(x$noise.sig)==1){
        	noise.sig<-x$noise.sig
		}else{
			warning("Length x$noise.sig > 1. Selecting first element.\n")
				noise.sig<-x$noise.sig[1]
		}
	}

	orig.ylim<-ylim
	main.orig<-main

	S0<-as.data.frame(x$S0,stringsAsFactors=FALSE)
	S0$Trial<-x$data[,trial.cn]

	if(length(trial)>1){
		warning("Length of trial > 1. Selecting first element.\n")
	}
	trial<-trial[1]

	if(length(ic)>1){
		warning("Length of ic > 1. Selecting first element.\n")
	}
	ic<-ic[1]

	# set-up temporary source matrix with nrow = one epoch
	S.temp<-S0[S0[,trial.cn]==trial,] 
	S.temp<-S.temp[,-which(colnames(S.temp)=="Trial")]
	S.temp<-as.matrix(S.temp)

	# set-up temporary x with nrow = one epoch
	x.temp=x$data[x$data[,trial.cn]==trial,]

	# get correlation
	cor.info<-x$correction.info
	cor.info<-cor.info[cor.info$NoiseSignal==noise.sig,]
	correlation.smry<-cor.info[cor.info$IC==ic&cor.info$Trial==trial,"Corr"]
	correlation<-cor(S.temp[,ic],x.temp[,noise.sig])
	if(abs(correlation)>=threshold){
		#if(round(correlation.smry,1)!=round(correlation,1))stop(paste("Problem with correlation between IC ",ic," and noise signal ",noise.sig," at trial ",trial,". Correlation from summary does not match computed correlation: Correlation summary = ",correlation.smry,"; Computed correlation = ",correlation,".\n",sep=""))
		if(length(col)==2){
			mycol<-col[2]
		}else{
			mycol<-col
		}
	}else{
		mycol<-"grey"
	}
	
	if(is.null(main.orig)){
		main<-paste(trial.cn," ",trial," -- Noise Channel ",noise.sig,
			" -- IC ",ic,sep="")
	}
	tmp<-scale(x.temp[,noise.sig])
	if(correlation<0){
		S.temp<-S.temp*-1
	}
	if(is.null(ylim)){
		ylim<-range(c(scale(S.temp[,ic]),tmp))
	}
	
	if(is.null(orig.ylim)){
		ylim<-range(c(tmp,S.temp[,ic]),na.rm=TRUE)
		ylim<-c(-1*max(abs(ylim)),max(abs(ylim)))
	}
	plot(tmp,type="l",col=col[1],main=main,xlab=xlab,
    	ylab=ylab,ylim=ylim,cex=cex,cex.main=cex,cex.lab=cex,
    	cex.axis=cex,...)
	lines(scale(S.temp[,ic]),col=mycol,cex=cex,...)
	# add correlation
    mtext(paste("correlation = ",round(correlation,3),sep=""),
    	side=3,line=0,adj=0,col=mycol,cex=cex*0.75,...)
	legend("topleft",legend=c("noise channel","IC"),
		col=c(col[1],mycol),bty="n",cex=cex,lty=1,...)
}

