Simplexplot <-
function(OBJ,items,main, ...){

	par(ask=TRUE)
	
	if(missing(main)){main=-1}

	plotit<-function(x,OBJ,main,...){

		Torank<-OBJ$OCC[which(OBJ$OCC[,1]==x),]
		
		order<-rank(apply(Torank[,-c(1:3)],1,max))
		order<-max(order)+1-order	
		
	
		highest3<-sapply(order,function(yyy)ifelse(yyy>3,0,1))
		
		subset<-which(highest3==max(highest3))
		
		Normed<-Torank[subset,-c(1:3)]

		ToPlot<-t(Normed)
		
		tonamecols<-as.character(Torank[subset,2])
		tonamecols[which(tonamecols==-1)]='NA'
		
		
			
		
		colnames(ToPlot)<-sapply(tonamecols,function(x)paste("    Option # :   ",x,sep=""))
		
		
		ToPlot<-ToPlot/apply(ToPlot,1,sum)
		
		pts<-length(OBJ$evalpoints)
		onethird<-ceiling(pts/3)

		if(main==-1){main=paste("Item : ",OBJ$itemlabels[[x]])}

		triax.plot(ToPlot,pch=c(rep(3,onethird),rep(0,onethird),rep(19,onethird)),col.symbols=c(rep("red",onethird),rep("green",onethird),rep("blue",onethird)),show.grid=TRUE,main=main,...)
		legend(x=1,y=1,c("Low Score", "Medium Score", "High Score"), pch=c(3,0,19),col=c("red","green","blue"))
	
		}
	
	nada<-sapply(items,plotit,OBJ=OBJ,main,...)
	

}

