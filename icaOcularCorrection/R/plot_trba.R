plot_trba<-function(x,data,channel,noise.sig=NULL,n.win=10,
				   new.page=TRUE,trial.cn="Trial",...){
	if(!"icac"%in%class(x))stop("object \"x\" is not of class \"icac\"\n")
	if(missing(x))stop("please supply an \"icac\" object \"x\"\n")
	if(missing(data))stop("Please provide an object \"data\" containing the original uncorrected data.\n")
	if(missing(channel))stop("please supply an channel\n")
	if(!trial.cn%in%colnames(data))stop(paste(trial.cn,"not in \"data\" column names\n"))

	cat("\n")
	cat("--------------------\n")
	cat("LEGEND:\n")
	cat("    * black line is uncorrected signal;\n")
	cat("    * blue line is corrected signal;\n")
	cat("    * grey line is noise signal.\n")
	cat("--------------------\n")
	cat("\n")
		
	# get trials
	trials<-sort(unique(data[,trial.cn]))

	odanp<-par()$ask
	devAskNewPage(new.page)

	# compute layout
	lm<-list()
	lm[[1]]<-c(1,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
	if(n.win>1){
		for(w in 2:n.win){
			lm[[w]]<-lm[[w-1]]+2
		}
	}
	lm<-unlist(lm)
	layout(matrix(lm,nrow=n.win,ncol=15,byrow=TRUE))	

	# plot
	for(tr in trials){
		par(mar=c(0.1,0.1,0.1,0.1))
		plot(0,0,ann=F,axes=F,type="n")
		text(x=0,y=0,labels=tr,bty="n",cex=1,col="red")

		od<-scale(data[data[,trial.cn]==tr,channel[1]])
		cd<-scale(x$data[x$data[,trial.cn]==tr,channel[1]])

		if(!is.null(noise.sig)){
			ns<-scale(data[data[,trial.cn]==tr,noise.sig[1]])
		}

		if(!is.null(noise.sig)){
			ylim<-c(od,cd,ns)
		}else{
			ylim<-c(od,cd)
		}
		ylim<-max(abs(ylim))
		ylim<-c(-1*ylim,ylim)

		plot(od,type="l",ann=FALSE,axes=FALSE,ylim=ylim,...)
		lines(cd,type="l",col="blue",...)

		if(!is.null(noise.sig)){
			lines(ns,type="l",col="gray",...)
		}
	}

	layout(matrix(1,nrow=1,ncol=1))
	par(mar=c(5.1,4.1,4.1,2.1))
	devAskNewPage(odanp)
}
