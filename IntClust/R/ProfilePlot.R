ProfilePlot<-function(Genes,Comps,GeneExpr=NULL,Raw=FALSE,OrderLab=NULL,ColorLab=NULL,nrclusters=NULL,cols=NULL,AddLegend=TRUE,margins=c(8.1,4.1,1.1,6.5),extra=5,plottype="new",location=NULL,...){
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	
	
	if(class(GeneExpr)[1]=="ExpressionSet"){
		GeneExpr <- Biobase::exprs(GeneExpr)
		
	}
	
	if(!is.null(OrderLab)){
		if(class(OrderLab)=="character"){
			orderlabs=OrderLab
		}
		else{
			orderlabs=OrderLab$Clust$order.lab
			GeneExpr=GeneExpr[,match(orderlabs,colnames(GeneExpr))]
		}
	}
	else{
		orderlabs=colnames(GeneExpr)
	}
	
	if(!is.null(ColorLab)){
		Data1 <- ColorLab$Clust
		ClustData1=stats::cutree(Data1,nrclusters) 
		
		ordercolors=ClustData1[Data1$order]
		names(ordercolors)=Data1$order.lab
		
		ClustData1=ClustData1[Data1$order]	
		
		order=seq(1,nrclusters)
		
		for (k in 1:length(unique(ClustData1))){
			select=which(ClustData1==unique(ClustData1)[k])
			ordercolors[select]=order[k]
		}
		
		if(!is.null(OrderLab)){
			if(class(OrderLab)=="character"){
				ordernames=OrderLab
			}
			else{
				ordernames=OrderLab$Clust$order.lab
			}	
			ordercolors=ordercolors[ordernames]
		}
		
		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)	
		
	}
	else{
		colors1<-rep("green",length(Comps))
		colors2<-rep("black",length(orderlabs[which(!(orderlabs%in%Comps))]))
		colors=c(colors1,colors2)
		AllCpds=c(Comps,orderlabs[which(!(orderlabs%in%Comps))])
		names(colors)=AllCpds
	}
	
	#yvalues=c()
	#allvalues=c()
	#for(i in 1:length(Genes)){
	#	yvalues=as.vector(GeneExpr[which(rownames(GeneExpr)==Genes[i]),])
	#	allvalues=c(allvalues,yvalues-mean(yvalues))
	#}	
	yvalues=GeneExpr[Genes,]
	if(Raw==FALSE & class(yvalues) != "numeric"){
		allvalues=as.vector(apply(yvalues,1,function(c) c-mean(c)))
	}
	else if(Raw==FALSE & class(yvalues) == "numeric"){
		allvalues=as.vector(sapply(yvalues,function(c) c-mean(yvalues)))
	}
	else{
		allvalues=as.vector(yvalues)
	}	
	ylims=c(min(allvalues)-0.1,max(allvalues)+0.1)
	
	plottypein(plottype,location)
	graphics::par(mar=margins,xpd=TRUE)
	graphics::plot(type="n",x=0,y=0,xlim=c(0,ncol(GeneExpr)),ylim=ylims,ylab=expression(log[2] ~ paste("fold ", "change")),xlab=" ",xaxt="n")
	#ylims=c()
	Indices=c(colnames(GeneExpr)[which(colnames(GeneExpr)%in%Comps)],colnames(GeneExpr)[which(!colnames(GeneExpr)%in%Comps)])
	
	for(i in 1:length(Genes)){
		GenesComps=as.numeric(GeneExpr[which(rownames(GeneExpr)==Genes[i]),colnames(GeneExpr)%in%Comps])
		Others=as.numeric(GeneExpr[which(rownames(GeneExpr)==Genes[i]),!(colnames(GeneExpr)%in%Comps)])
		if(length(Others)==0){
			Continue=FALSE
		}else{Continue=TRUE}
		
		#ylims=c(ylims,c(GenesComps,Others)-mean(c(GenesComps,Others)))	
		if(Raw==FALSE){
			yvalues1=GenesComps-mean(c(GenesComps,Others))	
		}
		else{
			yvalues1=GenesComps
		}
		
		graphics::lines(x=seq(1,length(GenesComps)),y=yvalues1,lty=1,col=i,lwd=1.6)
		#points(x=seq(1,length(GenesComps)),y=yvalues1,pch=19,col=i)
		graphics::segments(x0=1,y0=mean(yvalues1[1:length(GenesComps)]),x1=length(GenesComps),y1=mean(yvalues1[1:length(GenesComps)]),lwd=1.5,col=i)
		
		
		if(Continue==TRUE){
			if(Raw==FALSE){
				yvalues2=Others-mean(c(GenesComps,Others))	
			}
			else{
				yvalues2=Others
			}
			
			
			graphics::lines(x=seq(length(GenesComps)+1,ncol(GeneExpr)),y=yvalues2,lty=1,col=i,lwd=1.6)
			graphics::segments(x0=length(GenesComps)+1,y0=mean(yvalues2[1:length(Others)]),x1=ncol(GeneExpr),y1=mean(yvalues2[1:length(Others)]),lwd=1.5,col=i)
			
		}

		
	}	
	#Indices=c(colnames(GeneExpr)[which(colnames(GeneExpr)%in%Comps)],colnames(GeneExpr)[which(!colnames(GeneExpr)%in%Comps)])
	if(!is.null(ColorLab)){
		graphics::axis(1, labels=FALSE)
		#box("outer")
		graphics::mtext(substr(Indices,1,15), side = 1,  at=seq(0.5,(ncol(GeneExpr)-0.5)), line=0.2, las=2, cex=0.70,col=colors[Indices])
		#mtext(substr(Indices,1,15), side = 1,  at=seq(0.5,(ncol(GeneExpr)-0.5)), line=0.2, las=2, cex=0.70,col=c(rep("blue",7),rep("black",(56-7))))
	}
	else{
		#axis(1,at=seq(0.5,(ncol(GeneExpr)-0.5)),labels=Indices,las=2,cex.axis=0.70,xlab=" ",col=colors[Indices])
		#axis(1,at=seq(0.5,(ncol(GeneExpr)-0.5)), labels=FALSE)
		graphics::mtext(Indices, side = 1,  at=seq(0.5,(ncol(GeneExpr)-0.5)),line=0.2, las=2, cex=0.70,col=colors[Indices])
	}
	graphics::axis(2,ylab=expression(log[2] ~ paste("fold ", "change")))
	
	if(AddLegend==TRUE){
		
		labels=Genes
		colslegend=seq(1,length(Genes))
		
		graphics::par(xpd=T,mar=margins)
		graphics::legend(ncol(GeneExpr)+extra,mean(c(min(ylims),max(ylims))),legend=c(labels),col=c(colslegend),lty=1,lwd=3,cex=0.8)
		
	}
	plottypeout(plottype)
}
