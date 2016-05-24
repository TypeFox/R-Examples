OCCplot <-
function(OBJ,items,alpha,axis,quants,main,xlab,ylab,xlim,ylim,...){


		if(missing(ylim)){ylim=c(0,1)}
		if(missing(xlim)){xlim=c(min(axis),max(axis))}
		if(missing(ylab)){ylab="Probability"}
		if(missing(main)){main=-1}	
		if(missing(alpha)){alpha<-FALSE}
		

	plotit<-function(x,OBJ,alpha,axis,quants,main, xlim,ylim,xlab,ylab,...){

		if(main==-1){main=paste("Item: ",OBJ$itemlabels[x],"\n")}


		plot(1,ylim=ylim,main=main,xlim=xlim,type="l",ylab=ylab,xlab=xlab,...)

		IRFlines<-OBJ$OCC[which(OBJ$OCC[,1]==x),]
		SE<-OBJ$stderrs[which(OBJ$OCC[,1]==x),]

		for(i in 1:nrow(IRFlines)){
			if(OBJ$scale[x]==1){
				if(IRFlines[i,3]==1){colortouse<-"blue"; lwidth <-2}
				else{colortouse<-"red"; lwidth <-1}
			}
			else{
				colortouse<-"black"; lwidth <-1;
			}
			
			lines(axis,IRFlines[i,-c(1:3)],col=colortouse, lwd=lwidth)
			word<-ifelse(IRFlines[i,2]==-1,"NA",as.character(IRFlines[i,2]))
			wordloc<-round(runif(1,min=10,max=OBJ$nevalpoints-10))
			text(axis[wordloc],IRFlines[i,wordloc],word,cex=.7)

			## Cint
			if(alpha){

				ME<-qnorm(1-alpha/2)*SE[i,-c(1:3)]

				confhigh<-sapply(IRFlines[i,-c(1:3)]+ME,function(x)min(x,1));
				conflow<-sapply(IRFlines[i,-c(1:3)]-ME,function(x)max(x,0));
		
				lines(axis,confhigh,lty=2,col=colortouse)
				lines(axis,conflow,lty=2,col=colortouse)

			

			}		
	

		}




		

		

			axis(3,at=quants, lab=labels(quants),tck=0)
			abline(v=quants,col="blue",lty=2)
		
		box()
	

	}
	par(ask=TRUE)

	nada<-sapply(items,plotit,OBJ=OBJ,alpha=alpha,axis=axis,quants=quants,main,xlim,ylim,xlab,ylab,...)
	

}

