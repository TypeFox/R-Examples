ICCDIFplot <-
function(OBJ,items,alpha,axis,quants,main, xlab,ylab,xlim,ylim,cex,...){

	


		if(missing(ylim)){ylim=-1}
		if(missing(cex)){cex=.4}
		if(missing(xlim)){xlim=c(min(axis),max(axis))}
		if(missing(ylab)){ylab="Expected Item Score"}
		if(missing(main)){main=-1}	
		if(missing(alpha)){alpha<-FALSE}


	
	plotit <- function(x,OBJ,alpha,axis,quants,main, xlab,ylab,xlim,ylim,cex,...){


		maxitem<-max(OBJ$OCC[which(OBJ$OCC[,1]==x),3])
		minitem<-min(OBJ$OCC[which(OBJ$OCC[,1]==x),3])

		if(main==-1){main=paste("Item: ",OBJ$itemlabels[x],"\n")}
		if(ylim==-1){ylim=c(minitem,maxitem)}


			Estimate0<-OBJ$OCC[which(OBJ$OCC[,1]==x),]
			Estimate1<-apply(Estimate0[,-c(1:3)],2,function(yyy)yyy*Estimate0[,3])
			EstimateFULL<-apply(Estimate1,2,sum)


			

		plot(c(min(axis),max(axis)),c(0,max(OBJ$OCC[which(OBJ$OCC[,1]==x),])),type="n",xlim=xlim, ylim=ylim ,xlab=xlab,ylab=ylab,main = main,...)
		
	

		
		lines(axis,EstimateFULL,col="black",lwd=1.5,...)

		ngrps<-length(OBJ$groups)
		plot_colors <- c("blue","red","forestgreen","black","yellow","orange")
		if(ngrps>6){plot_colors<-c(plot_colors,sample(colors(),ngrps-6))}
		
	


		for(i in 1:ngrps){

			cgrp<-OBJ$DIF[[i]]



			dbins<-cut(cgrp$subjtheta,breaks=c(-999,cgrp$evalpoints[-length(cgrp$evalpoints)],999),labels=FALSE)

			resp0<-cgrp$binaryresp[which(cgrp$binaryresp[,1]==x),]
			respit1<-apply(resp0[,-c(1:3)],2,function(x)x*resp0[,3])
			respit<-apply(respit1,2,sum)
			#respit<-resp0[which(resp0[,3]==maxitem),-c(1:3)]
			propevalpoints<-numeric()

			for (jjj in 1:cgrp$nevalpoints){
		
				binaryrespp<-respit[which(dbins==jjj)]
				propevalpoints[jjj]<-sum(binaryrespp)/length(binaryrespp)
			}

			Estimate0<-cgrp$OCC[which(cgrp$OCC[,1]==x),]
			Estimate1<-apply(Estimate0[,-c(1:3)],2,function(yyy)yyy*Estimate0[,3])
			Estimate<-apply(Estimate1,2,sum)





			Stderr0<-cgrp$stderrs[which(cgrp$OCC[,1]==x),]
			Stderr1<-apply(Stderr0[,-c(1:3)],2,function(zzz)zzz*Stderr0[,3])
			Stderr<-apply(Stderr1,2,sum)
			
		
			lines(axis,Estimate,col=plot_colors[i], ...)
			points(axis,propevalpoints,cex=cex,col=plot_colors[i], ...)



			if(alpha){		
				SE<-try(qnorm(1-alpha/2)*Stderr)
	

				confhigh<-try(sapply(Estimate+SE,function(x)min(x,maxitem)))
				conflow<-try(sapply(Estimate-SE,function(x)max(x,minitem)))


				try(lines(axis,confhigh,lty=2,col=plot_colors[i]))
				try(lines(axis,conflow,lty=2,col=plot_colors[i]))
			}



		}

	

		#resp0<-OBJ$responses[which(OBJ$responses[,1]==x),]
		#respit<-resp0[which(resp0[,3]==maxitem),-c(1:3)]

	
	axis(3,at=quants, lab=labels(quants),tck=0)
	legend(min(axis), maxitem,c("Overall",OBJ$groups), cex=0.8, col=c("black",plot_colors[1:ngrps]),lwd=2, bty="n");
	abline(v=quants,col="blue",lty=2)
	box()
	
	}

	
	par(ask=TRUE)

	nada<-sapply(items,plotit,OBJ=OBJ,alpha=alpha,axis=axis,quants=quants,main,xlab,ylab,xlim,ylim,cex,...)
	

}

