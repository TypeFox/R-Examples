ICCplot <-
ICCplot <-
function(OBJ,items,alpha,axis,quants,main, xlab,ylab,xlim,ylim,cex,...){

	dbins<-cut(OBJ$subjtheta,breaks=c(-999,OBJ$evalpoints[-length(OBJ$evalpoints)],999),labels=FALSE)


		if(missing(ylim)){ylim=-1}
		if(missing(cex)){cex=.4}
		if(missing(xlim)){xlim=c(min(OBJ$evalpoints),max(OBJ$evalpoints))}
		if(missing(ylab)){ylab="Expected Item Score"}
		if(missing(main)){main=-1}	
		if(missing(alpha)){alpha<-.05}



	plotit <- function(x,OBJ,alpha,axis,quants,main, xlab,ylab,xlim,ylim,cex,...){


	Estimate0<-OBJ$OCC[which(OBJ$OCC[,1]==x),]
	maxitem<-max(Estimate0[,3])
	minitem<-min(Estimate0[,3])
	Stderr0<-OBJ$stderrs[which(OBJ$OCC[,1]==x),]
	resp0<-OBJ$binaryresp[which(OBJ$binaryresp[,1]==x),]


#	Estimate<-Estimate0[which(Estimate0[,3]==maxitem),]
#	Stderr<-Stderr0[which(Estimate0[,3]==maxitem),]

	Estimate1<-apply(Estimate0[,-c(1:3)],2,function(x)x*Estimate0[,3])
	Estimate<-apply(Estimate1,2,sum)
	#Estimate <- sapply(Estimate,function(x) min(0,max(1,x)));

	Stderr1<-apply(Stderr0[,-c(1:3)],2,function(x)x*Stderr0[,3])
	Stderr<-apply(Stderr1,2,sum)



	respit1<-apply(resp0[,-c(1:3)],2,function(x)x*resp0[,3])
	respit<-apply(respit1,2,sum)
	
#	respit<-resp0[which(resp0[,3]==maxitem),-c(1:3)]
	propevalpoints<-numeric()

	for (i in 1:OBJ$nevalpoints){
		
		binaryrespp<-respit[which(dbins==i)]
		propevalpoints[i]<-sum(binaryrespp)/length(binaryrespp)
	}


		if(main==-1){main=paste("Item: ",OBJ$itemlabels[x],"\n")}
		if(ylim==-1){ylim=c(minitem,maxitem)}

		plot(axis,Estimate,ylim=ylim,type="l",ylab=ylab,xlab=xlab,main=main,...)
		
		 
		 
		
		
		if(alpha){

			SE<-qnorm(1-alpha/2)*Stderr

			
			confhigh<-sapply(Estimate+SE,function(x)min(x,maxitem));
			conflow<-sapply(Estimate-SE,function(x)max(x,minitem));
			#print(Estimate-SE)
			#print(conflow)
			
		

			lines(axis,confhigh,lty=2,col="red")
			lines(axis,conflow,lty=2,col="red")

		}

		points(axis,propevalpoints,cex=cex,...)	
		
		
			axis(3,at=quants, lab=labels(quants),tck=0)
			abline(v=quants,col="blue",lty=2)
		

		box()

		return(Estimate)
		
	}

	
	par(ask=TRUE)

	sapply(items,plotit,OBJ=OBJ,alpha=alpha,axis=axis,quants=quants,main,xlab,ylab,xlim,ylim,cex,...)
	

}

