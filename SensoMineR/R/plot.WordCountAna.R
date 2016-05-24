plot.WordCountAna<-function(x,axes=c(1,2),choix="prod",lab=TRUE,color=NULL,pch=NULL,proba=0.05,xlim=NULL,ylim=NULL,cex=1,title=NULL,new.plot=TRUE,...){
	res.WCA<-x
	if (!inherits(res.WCA,"WordCountAna")) stop("non convenient data")
	if (choix=="prod"){
		if (is.null(title)) title<-"Products"
        	coord<-res.WCA$mfact$ind$coord[,axes,drop=FALSE]
	  	if (is.null(color)) color="blue"
	  	if (is.null(pch)) pch=16
	}
	if (choix=="panel"){
		if (is.null(title)) title<-"Panellists"
        	coord<-res.WCA$mfact$group$coord[,axes,drop=FALSE]
        	if (is.null(color)) color="brown"
	  	if (is.null(pch)) pch=17
	  	xlim<-ylim<-c(0,1)
	}
	if (choix=="dist"){
	  if (is.null(title)) title<-"Distinct-words"
        coord<-res.WCA$centroids[,axes,drop=FALSE]
        if (is.null(color)) color="black"
	  if (is.null(pch)) pch=18
	}
	if (choix=="cons"){	
		if (is.null(title)) title<-"Consensual words"
        	coord<-res.WCA$centroids[which(rownames(res.WCA$centroids)%in%rownames(res.WCA$cons)[res.WCA$cons[,3]<proba]),axes,drop=FALSE]
        	if (is.null(color)) color="red"
	  	if (is.null(pch)) pch=15
	}
	lab.x<-paste("Dim ",axes[1]," (",signif(res.WCA$mfact$eig[axes[1],2],4)," %)",sep="")
	lab.y<-paste("Dim ",axes[2]," (",signif(res.WCA$mfact$eig[axes[2],2],4)," %)",sep="")
    	if (new.plot) dev.new()
	if (is.null(xlim)) {
            xmin<-xmax<-0
            xmin<-min(xmin,coord[,1])
            xmax<-max(xmax,coord[,1])
		xlim<-c(xmin,xmax)*1.2
	}
	if (is.null(ylim)) {
   		ymin<-ymax<-0
         	ymin<-min(ymin,coord[,2])
         	ymax<-max(ymax,coord[,2])	
		ylim<-c(ymin,ymax)*1.2
	}
	plot(coord,xlab=lab.x,ylab=lab.y,xlim=xlim,ylim=ylim,pch=pch,col=color,cex=cex,main=title,asp=1,...)
	if (lab) text(coord[,1],coord[,2],labels=rownames(coord),pos=3,cex=cex,col=color)
	if (choix!="panel") abline(v=0,h=0,lty=2)
}