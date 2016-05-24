BoxPlotDistance<-function(Data1,Data2,type=c('data','dist','clusters'),distmeasure="tanimoto",normalize=FALSE,method=NULL,lab1,lab2,limits1=NULL,limits2=NULL,plot=1,StopRange=FALSE,plottype="new",location=NULL){
	type<-match.arg(type)
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
	
	C1=D1=C2=D2=NULL
	if(type=='clusters'){
	
		Dist1<-Data1$DistM
		Dist2<-Data2$DistM
	

	}
	else if(type=='data'){
		Dist1<-Distance(Data1,distmeasure[1],normalize,method)
		Dist2<-Distance(Data2,distmeasure[2],normalize,method)
		DistL=list(Dist1,Dist2)
		for(i in 1:2){
			if(StopRange==FALSE & !(0<=min(DistL[[i]]) & max(DistL[[i]])<=1)){
				message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				DistL[[i]]=Normalization(DistL[[i]],method="Range")
			}
		}
		
	}
	else if(type=='dist'){
		Dist1=Data1
		Dist2=Data2	
		DistL=list(Dist1,Dist2)
		for(i in 1:2){
			if(StopRange==FALSE &  !(0<=min(DistL[[i]]) & max(DistL[[i]])<=1)){
				message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				DistL[[i]]=Normalization(DistL[[i]],method="Range")
			}
		}
	}
	
	OrderNames=rownames(Dist1)
	Dist2=Dist2[OrderNames,OrderNames]
	
	Dist1lower <- Dist1[lower.tri(Dist1)]
	Dist2lower <- Dist2[lower.tri(Dist2)]
	
	Categorize<-function(Distlower,limits){
		Cat=c(rep(0,length(Distlower)))
		for(j in 1:(length(limits)+1)){
			if(j==1){
				Cat[Distlower<=limits[j]]=j
			}
			else if(j<=length(limits)){
				Cat[Distlower>limits[j-1] & Distlower<=limits[j]]=j
			}	
			else{
				Cat[Distlower>limits[j-1]]=j
			}
		}
		Cat<-factor(Cat)
		return(Cat)
		
	}
	
	#plot2
	if(!(is.null(limits1))){
		Dist1cat<-Categorize(Dist1lower,limits1)
		
		dataBox2<-data.frame(D2=Dist2lower,C1=Dist1cat)
		p2<-ggplot2::ggplot(dataBox2,ggplot2::aes(factor(C1),D2)) #x,y
		p2<-p2+ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_point(color="blue",size=2,shape=19,position="jitter",cex=1.5)+
				ggplot2::xlab(lab1)+ggplot2::ylab(lab2)
	}
	#plot1
	if(!(is.null(limits2))){
		Dist2cat<-Categorize(Dist2lower,limits2)
		dataBox1<-data.frame(D1=Dist1lower,C2=Dist2cat)
		p1<-ggplot2::ggplot(dataBox1,ggplot2::aes(factor(C2),D1)) #x,y
		p1<-p1+ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_point(color="blue",size=2,shape=19,position="jitter",cex=1.5)+
				ggplot2::xlab(lab2)+ggplot2::ylab(lab1)
		
	}	
	if(type==3){
		if(plottype=="pdf"){
			location=paste(location,'_type3.pdf',sep="")
		}
		plottypein(plottype,location)
		gridExtra::grid.arrange(p1, p2, ncol=2,nrow=1)		

	}
	else if(type==1){
		if(plottype=="pdf"){
			location=paste(location,'_type1.pdf',sep="")
		}
		plottypein(plottype,location)
		print(p1)
		
	}
	else if(type==2){
		if(plottype=="pdf"){
			location=paste(location,'_type2.pdf',sep="")
		}
		plottypein(plottype,location)
		print(p2)		
	}
	
}
