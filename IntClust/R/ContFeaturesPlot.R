ContFeaturesPlot<-function(LeadCpds,Data,nrclusters=NULL,OrderLab=NULL,ColorLab=NULL,cols=NULL,ylab="bio-assays",AddLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),plottype="new",location=NULL){
	
	if(all(LeadCpds%in%rownames(Data))){
		Data=t(Data)
	}
	
	if(!is.null(OrderLab)){
		if(class(OrderLab)=="character"){
			orderlabs=OrderLab
		}
		else{
			orderlabs=OrderLab$Clust$order.lab
			Data=Data[,match(orderlabs,colnames(Data))]
		}
	}
	else{
		orderlabs=colnames(Data)
	}
	
	temp=orderlabs[which(!(orderlabs%in%LeadCpds))]
	AllCpds=c(LeadCpds,temp)
	Data=Data[,AllCpds]
	
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
	

	plottypein(plottype,location)
	graphics::par(mar=margins)
	graphics::plot(x=0,y=0,xlim=c(0,(ncol(Data)+3)),ylim=c(min(Data)-0.5,max(Data)+0.5),type="n",ylab=ylab,xlab='',xaxt='n')
	for(i in c(1:nrow(Data))){	
		graphics::lines(x=seq(1,length(LeadCpds)),y=Data[i,which(colnames(Data)%in%LeadCpds)],col=i)
	}
	for(i in c(1:nrow(Data))){	
		graphics::lines(x=seq(length(LeadCpds)+4,(ncol(Data)+3)),y=Data[i,which(!(colnames(Data)%in%LeadCpds))],col=i)
	}

	

	if(!(is.null(ColorLab)) | is.null(nrclusters)){
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
	
		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)	
	}
	else{
		colors1<-rep("green",length(LeadCpds))
		colors2<-rep("black",length(temp))
		colors=c(colors1,colors2)
		names(colors)=AllCpds
	}
	
	graphics::mtext(LeadCpds,side=1,at=seq(1,length(LeadCpds)),line=0,las=2,cex=0.70,col=colors[LeadCpds])
	graphics::mtext(temp, side = 1, at=c(seq(length(LeadCpds)+4,(ncol(Data)+3))), line=0, las=2, cex=0.70,col=colors[temp])
	if(AddLegend==TRUE){
		
		labels=rownames(Data)
		colslegend=seq(1,length(rownames(Data)))
		
		graphics::par(xpd=T,mar=margins)
		graphics::legend(ncol(Data)+5,mean(c(min(c(min(Data)-0.5,max(Data)+0.5)),max(c(min(Data)-0.5,max(Data)+0.5)))),legend=c(labels),col=c(colslegend),lty=1,lwd=3,cex=0.8)
		
	}
	plottypeout(plottype)
}