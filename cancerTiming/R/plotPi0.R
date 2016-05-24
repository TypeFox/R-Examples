plotPi0<-function(est,ui,li,samples,labels=NULL,xlab="Segment",ylab=expression(pi[0]),pchFac=rep(1,length(est)), 
	whHighlight=NULL,nMut=NULL,xorder=NULL,...){
	#require(gplots)
	if(requireNamespace("gplots", quietly = TRUE) ){
		samples<-factor(samples)
		indToPlot<-tapply(1:length(est),samples, function(ii){
			p<-est[ii]
			whNA<-ii[which(is.na(p))]
			return(setdiff(ii[order(p)],c(whNA)))
			#whL20<-ii[which(df$N[ii]<20)]
			#return(setdiff(ii[order(p)],c(whNA,whL20)))
		})
		if(is.null(xorder)) xord<-unlist(indToPlot)
		else xord<-xorder
		nPerSamp<-sapply(indToPlot,length)
		plotCI(est[xord],ui=ui[xord],li=li[xord],xlab=xlab,
			ylab=ylab,xaxt="n",pch=c(1:5)[pchFac[xord]],sfrac=0.001,cex=.7,...)
		n<-length(xord)
	
		#separate the samples
		segments(x0=cumsum(nPerSamp)+.5,y1=par("usr")[4]+.2,y0=par("usr")[3]-.1,xpd=NA,lty=2)
		mids<-(c(1,head(cumsum(nPerSamp)+.5,-1))+cumsum(nPerSamp)+.5)/2
		axis(3,mids,levels(samples),tick=FALSE)
	
		#label the segments
		axis(1,1:n,labels[xord],las=2,cex.axis=0.7)	
		if(!is.null(nMut)){
			axis(3,1:n,nMut[xord],line=-1,tick=FALSE,cex.axis=.7)	
			 mtext("N=",side=3,at=par("usr")[1])		
		}
	
		#
		if(!is.null(whHighlight)){
			wh<-which(xord%in%whHighlight)
			plotCI((1:n)[wh],est[xord][wh],ui=ui[xord][wh],li=li[xord][wh],
				pch=c(1:5)[pchFac[xord][wh]],col="red",sfrac=0.001,cex=.7,add=TRUE,lwd=1.7)
		
		}
		invisible(xord)
	}
	else{stop("This function requires the package 'gplots' to be installed.")}
}