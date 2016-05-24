OCCDIFplot <-function(OBJ,items,alpha,axis,quants,main,xlab,ylab,xlim,ylim,...){

		if(missing(ylim)){ylim=c(0,1)}
		if(missing(xlim)){xlim=c(min(axis),max(axis))}
		if(missing(ylab)){ylab="Probability"}
		if(missing(main)){main0=-1}	
		if(missing(alpha)){alpha<-FALSE}

	plotit<-function(x,OBJ,alpha,axis,quants,main, xlim,ylim,xlab,ylab,...){


		ngrps<-length(OBJ$groups)
		plot_colors <- c("blue","red","forestgreen","black","yellow","orange")
		if(ngrps>6){plot_colors<-c(plot_colors,sample(colors(),ngrps-6))}
		
		IRFlines<-OBJ$OCC[which(OBJ$OCC[,1]==x),]
		opts<-unique(IRFlines[,2])
		nopts<-length(opts)

		for(i in 1:nopts){
			
			if(main0==-1){main=paste("Item: ",OBJ$itemlabels[x]," ","Option: ",ifelse(opts[i]==-1,"NA",opts[i]),"\n\n")}
		


				plot(1,ylim=ylim,main=main,xlim=xlim,type="l",ylab=ylab,xlab=xlab,...)
	
			
			
			
			IRFlineBIG<-OBJ$OCC[which(OBJ$OCC[,1]==x & OBJ$OCC[,2]==opts[i]),]
			lines(axis,IRFlineBIG[-c(1:3)],col="black",lwd=1.5,...)
			
			for(j in 1:ngrps){

				cgrp<-OBJ$DIF[[j]]


				if(length(which(cgrp$OCC[,1]==x & cgrp$OCC[,2]==opts[i]))==0){next;}

				IRFline<-cgrp$OCC[which(cgrp$OCC[,1]==x & cgrp$OCC[,2]==opts[i]),]
				SE<-cgrp$stderrs[which(cgrp$OCC[,1]==x & cgrp$OCC[,2]==opts[i]),]
					
				try(lines(axis,IRFline[-c(1:3)],col=plot_colors[j],...))
				
			
				
	

			
				if(alpha){
					ME<-qnorm(1-alpha/2)*SE[-c(1:3)]
	
					confhigh<-sapply(IRFline[-c(1:3)]+ME,function(x)min(x,1))
					conflow<-sapply(IRFline[-c(1:3)]-ME,function(x)max(x,0))
		
					try(lines(axis,confhigh,lty=2,col=plot_colors[j]))
					try(lines(axis,conflow,lty=2,col=plot_colors[j]))

				}		
	

			}

			axis(3,at=quants, lab=labels(quants),tck=0)
				legend(min(axis), 1,c("Overall",OBJ$groups), cex=0.8, col=c("black",plot_colors[1:ngrps]),lwd=2, bty="n");
				abline(v=quants,col="blue",lty=2)
				box()

		}
		


	}

	par(ask=TRUE)

	nada<-sapply(items,plotit,OBJ=OBJ,alpha=alpha,axis=axis,quants=quants, main,xlim,ylim,xlab,ylab,...)
	

}

