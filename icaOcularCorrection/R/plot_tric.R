plot_tric<-function(x,ic,noise.sig=NULL,S="S0",trial.cn="Trial",trials=NULL,n.win=10,new.page=TRUE,...){
	if(!"icac"%in%class(x))stop("object not of class \"icac\"\n")
	if(missing(x))stop("please supply an \"icac\" object\n")
	#if(missing(noise.sig))stop("please supply a noise signal\n")
	if(missing(ic))stop("please supply an independent component\n")
	
	cat("\n")
	cat("--------------------\n")
	cat("LEGEND:\n")
	cat("    * black line is independent component;\n")
	cat("    * blue line is noise signal above threshold;\n")
	cat("    * grey line is noise signal below threshold.\n")
	cat("--------------------\n")
	cat("\n")

	threshold<-x$threshold
		
	if(is.null(trials)){
		trials<-sort(unique(x$data[,trial.cn]))
	}

	S<-as.data.frame(x[[S]])
	colnames(S)<-paste("IC",1:x$n.comp,sep="")
	S[,trial.cn]<-x$data[,trial.cn]
	
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
		tmp<-S[S[,trial.cn]==tr,paste("IC",ic,sep="")]
		if(length(unique(tmp))!=1){
			tmp<-scale(tmp)
		}

		if(!is.null(noise.sig)){
			tmp2<-scale(x$data[x$data[,trial.cn]==tr,noise.sig])
			if(cor(tmp,tmp2)<0){
				tmp2<-tmp2*-1
			}
		}

		if(exists("tmp2")){
			ylim<-range(c(tmp,tmp2))
		}else{
			ylim<-range(c(tmp))
		}
		
		plot(tmp,type="l",ann=FALSE,axes=FALSE,ylim=ylim)
		if(!is.null(noise.sig)){
			if(noise.sig%in%x$noise.sig){
				correlation<-cor(tmp,tmp2)
				if(abs(correlation)>=threshold){
					mycol<-"blue"
					if(correlation<0){
						tmp2<--1*tmp2
					}
				}else{
					mycol<-"grey"
					if(correlation<0){
						tmp2<--1*tmp2
					}
				}
				lines(tmp2,type="l",col=mycol)
			}else{
				stop("wrong noise.sig name.\n")
			}
		}
	}

	layout(matrix(1,nrow=1,ncol=1))
	par(mar=c(5.1,4.1,4.1,2.1))
	devAskNewPage(odanp)
}
