BinFeaturesPlot<-function(LeadCpds,OrderLab,Features,Data,ColorLab,nrclusters=NULL,cols=NULL,name=c("FP"),colors1=c('gray90','blue'),colors2=c('gray90','green'),margins=c(5.5,3.5,0.5,5.5),plottype="new",location=NULL){
	
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
	
	x<-c(1:length(AllCpds)) #x=comps
	y<-c(1:length(Features)) #y=feat
	PlotData<-t(Data[as.character(Features),AllCpds])
	plottypein(plottype,location)
	graphics::par(mar=margins)
	graphics::image(x,y,PlotData,col=colors1,xlab="",axes=FALSE,ann=FALSE,xaxt='n')
	if(length(unique(as.vector(PlotData[1:length(LeadCpds),])))==1){
		colors2=c("green")
	}
	graphics::image(x[1:length(LeadCpds)],y,PlotData[1:length(LeadCpds),],col=colors2,add=TRUE,xlab="",axes=FALSE,ann=FALSE,xaxt='n')
	
	if(!(is.null(ColorLab)) & !is.null(nrclusters)){
		if(attributes(ColorLab)$method=="Ensemble"){
			Data1<-ColorLab$Clust
			ClustData1=ColorLab$Clust$Clusters
		
		}
		else{
			Data1 <- ColorLab$Clust
			ClustData1=stats::cutree(Data1,nrclusters) 
		}
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
		colors1lab<-rep("green",length(LeadCpds))
		colors2lab<-rep("blue",length(temp))
		colors=c(colors1lab,colors2lab)
		names(colors)=AllCpds
	}
	
	graphics::mtext(colnames(PlotData), side = 4, at= c(1:ncol(PlotData)), line=0.2, las=2,cex=0.8)
	graphics::mtext(name, side = 2,  line=1, las=0, cex=1)
	graphics::mtext(rownames(PlotData), side = 1, at= c(1:nrow(PlotData)), line=0.2, las=2, cex=0.8,col=colors[AllCpds])
	plottypeout(plottype)
}
