plot_avgba<-function(x,data,channel=NULL,n.win=NULL,
				 new.page=TRUE,time.cn="Time",...){
	if(!"icac"%in%class(x))stop("object not of class \"icac\"\n")
	if(missing(x))stop("please supply an \"icac\" object\n")
	if(missing(data))stop("please supply the pre-correction data frame\n")

	cat("\n")
	cat("--------------------\n")
	cat("LEGEND:\n")	
	cat("    * black line is uncorrected signal;\n")
	cat("    * blue line is corrected signal.\n")
	cat("--------------------\n")	
	cat("\n")

	odanp<-par()$ask
	devAskNewPage(new.page)

	if(is.null(channel)){
		channel<-x$channel
	}

	# determine layout
	if(is.null(n.win)){
		if(length(channel)>10){
			n.win<-10
		}else{
			n.win=length(channel)
		}
	}
	lm<-list()
	lm[[1]]<-c(1,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
	if(length(channel)>1){
		for(w in 2:n.win){
			lm[[w]]<-lm[[w-1]]+2
		}
	}
	lm<-unlist(lm)
	layout(matrix(lm,nrow=n.win,ncol=15,byrow=TRUE))

	# plot
	for(e in channel){
		par(mar=c(0.1,0.1,0.1,0.1))
		plot(0,0,ann=F,axes=F,type="n")
		text(x=0,y=0,labels=e,bty="n",cex=1,col="red")

		od <- tapply(data[,e], data[,time.cn], mean)
		cd <- tapply(x$data[,e], x$data[,time.cn], mean)
		plot(as.numeric(names(od)),scale(od),type="l",ann=FALSE,axes=FALSE,...)
		lines(as.numeric(names(cd)),scale(cd),col="blue",...)
	}

	layout(matrix(1,nrow=1,ncol=1))
	par(mar=c(5.1,4.1,4.1,2.1))
	devAskNewPage(odanp)
}
