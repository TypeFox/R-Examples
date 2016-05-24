SelectnrClusters<-function(List,type=c("data","dist","pam"),distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,nrclusters = seq(5, 25, 1),names=NULL,StopRange=FALSE,plottype="new",location=NULL){
	
	type=match.arg(type)
	avsilwidth<-matrix(0,ncol=length(List),nrow=length(nrclusters))
	pamfunction<-function(DistM,nrclusters){
		asw=sapply(nrclusters, function(x) cluster::pam(DistM,x)$silinfo$avg.width)
		return(asw)
	}
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize,method))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		avsilwidth=sapply(Dist,function(x) pamfunction(x,nrclusters=nrclusters))
		rownames(avsilwidth)=nrclusters
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		avsilwidth=sapply(Dist,function(x) pamfunction(x,stats::start,stats::end))
		rownames(avsilwidth)=nrclusters
	}
	else{
		avsilwidth=sapply(List,function(x) return(x$silinfo$avg.width))
	}
	
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
	
	
	if(is.null(names)){
		names1=c()
		names2=c()
		for(i in 1:length(List)){
			names1=c(names1,paste("Silhouette widths for Data",i,sep=" "))
			names2=c(names2,paste("Nr Clusters for Data",i,sep=' '))
		}
		names1=c(names1,"Average Silhoutte Widths")
		names2=c(names2,"Optimal nr of clusters")
		
	}
	else{
		names1=c()
		names2=c()
		for(i in 1:length(List)){
			names1=c(names1,paste("Silhouette widths for",names[i],sep=" "))
			names2=c(names2,paste("Nr Clusters for",names[i],sep=' '))
		}
		names1=c(names1,"Average Silhoutte Widths")
		names2=c(names2,"Optimal nr of clusters")
	}
	

	rownames(avsilwidth)=nrclusters
	
	
	avsil=apply(avsilwidth,1,mean)
	avsilwidth=cbind(avsilwidth,avsil)
	colnames(avsilwidth)=names1
	
	
	plotsil<-function(sils,plottype,location,name){
		k.best=as.numeric(names(sils)[which.max(sils)])
		cat("silhouette-optimal number of clusters:", k.best, "\n")
		plottypein(plottype,location)
		graphics::plot(nrclusters, sils, type= "h", main = name,
				xlab= "k  (# clusters)", ylab = "average silhouette width")
		graphics::axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")
		plottypeout(plottype)
	}
	
	sapply(c(1:ncol(avsilwidth)),function(x) plotsil(avsilwidth[,x],plottype,location,names1[x]))
	
	Output=list()
	Output[[1]]=avsilwidth
	nrclusters=apply(avsilwidth,2,function(x) return(as.numeric(names(x)[which.max(x)])))
	nrclusters=as.data.frame(t(nrclusters))
	colnames(nrclusters)=names2
	rownames(nrclusters)="NrClusters"
	
	
	Output[[2]]=nrclusters
	
	names(Output)=c("Silhoutte_Widths","Optimal_Nr_of_CLusters")
	return(Output)
}
