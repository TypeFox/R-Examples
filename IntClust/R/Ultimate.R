Ultimate<-function(List,type=c("data","dist","clusters"),distmeasure,normalize=FALSE,method=NULL,StopRange=FALSE,NN=20,mu=0.5,T=20,t=10,r=NULL,nrclusters=NULL,nrclusterssep=c(7,7),nrclustersseq=NULL,weight=NULL,Clustweight=0.5,clust="agnes",linkage=c("ward","ward"),alpha=0.625,gap=FALSE,maxK=50,IntClust=c("ADC","ADECa","ADECb","ADECc","WonM","CECa","CECb","CECc","WeightedClust","WeightedSim","SNFa","SNFb","SNFc"),fusionsLog=TRUE,WeightClust=TRUE,PlotCompare=FALSE,cols=NULL,...){
	#Checking of conditions:
	type=match.arg(type)
	
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	

	
	#Separate hierarchical clustering on each data matrix
	
	SeparateClustering=list()
	
	for (i in 1:length(List)){
		SeparateClustering[[i]]=Cluster(List[[i]],type,distmeasure[i],normalize,method,clust,linkage,alpha,gap,maxK,StopRange)
		names(SeparateClustering)[i]=paste("Data",i,sep=" ")
	}
	
	
	#If FollowSel is True, cut tree into TibShirani number of cluster, pick last cluster number or select number of 
	#clusters and cluster choice yourself (recommended to perform separare clustering beforehand to make a selection)
	#nClust overrules optimal number of clusters determined by gap.
	
	if(is.null(nrclusters) & is.null(nrclusterssep) & gap==FALSE){
		stop("Please specify a number of clusters in nrclusters and optional in nrclusterssep")
	}
	
	if(!(is.null(nrclusters)) & is.null(nrclusterssep)){
		nrclusterssep=rep(nrclusters, length(List))		
	}
	
	if(is.null(nrclusters) & is.null(nrclusterssep) & gap==TRUE){
		nrclusterssep=c()
		for (i in 1:length(SeparateClustering)){
			temp=SeparateClustering[[i]]
			nrclusterssep[i]=temp$k[which(colnames(temp$k)=="Tibs2001SEmax")]
		}
		nrclusters=ceiling(mean(nrclusterssep))
	}
	
	if(is.null(nrclustersseq) & ( "ADECb" %in% IntClust | "WonM" %in% IntClust | "CECb" %in% IntClust) ){
		stop("Method ADECb, WonM and CECb need specification of nrclustersseq.")	
	}
	
	
	#Perform Integrated Cluster analysis of choice
	#SNF
	
	IntClustering=list()
	FollowedSelect=list()
	for (i in 1:length(IntClust)){
		if(IntClust[i]=="SNFa"){
			IntClustering[[i]]=SNFa(List,type,distmeasure,normalize,method,NN,mu,T,clust,linkage,alpha,StopRange)				
			names(IntClustering)[i]="SNFa"
			message("SNFa")
		}
		
		if(IntClust[i]=="SNFb"){
			IntClustering[[i]]=SNFb(List,type,distmeasure,normalize,method,NN,mu,T,clust,linkage,alpha,StopRange)				
			names(IntClustering)[i]="SNFb"
			message("SNFb")
		}
		
		if(IntClust[i]=="SNFc"){
			IntClustering[[i]]=SNFc(List,type,distmeasure,normalize,method,NN,mu,T,clust,linkage,alpha,StopRange)				
			names(IntClustering)[i]="SNFc"
			message("SNFc")
		}
		
		if(type=="data"){
		if(IntClust[i]=="ADC"){
			if(length(unique(distmeasure)) != 1){
				stop("When using ADClust, distance measures must be the same for all data matrices ")				
			}
			IntClustering[[i]]=ADC(List,unique(distmeasure),normalize,method,clust,linkage,alpha)
			names(IntClustering)[i]="ADC"
			message("ADC")
		}
		if(IntClust[i]=="ADECa"){
			if(length(unique(distmeasure)) != 1){
				stop("When using ADECa, distance measures must be the same for all data matrices ")				
			}
			IntClustering[[i]]=ADECa(List,unique(distmeasure),normalize,method,t,r,nrclusters,clust,linkage,alpha)
			names(IntClustering)[i]="ADECa"
			message("ADECa")
		}
		
		if(IntClust[i]=="ADECb"){
			if(length(unique(distmeasure)) != 1){
				stop("When using ADECb, distance measures must be the same for all data matrices ")				
			}
			IntClustering[[i]]=ADECb(List,unique(distmeasure),normalize,method,nrclusters=nrclustersseq,clust,linkage,alpha)
			names(IntClustering)[i]="ADECb"
			message("ADECb")
		}
		
		if(IntClust[i]=="ADECc"){
			if(length(unique(distmeasure)) != 1){
				stop("When using ADECc, distance measures must be the same for all data matrices ")				
			}
			IntClustering[[i]]=ADECc(List,unique(distmeasure),normalize,method,t,r,nrclusters=nrclustersseq,clust,linkage,alpha)
			names(IntClustering)[i]="ADECc"
			message("ADECc")
		}
		if(IntClust[i]=="CECa"){
			IntClustering[[i]]=CECa(List,distmeasure,normalize,method,t,r,nrclusters=nrclusterssep,weight,clust,linkage,alpha,Clustweight,StopRange)
			names(IntClustering)[i]="CECa"
			message("CECa")
		}
		
		if(IntClust[i]=="CECb"){
			
			IntClustering[[i]]=CECb(List,distmeasure,normalize,method,nrclusters=nrclustersseq,weight,clust,linkage,alpha,Clustweight,StopRange)
			names(IntClustering)[i]="CECb"
			message("CECb")
		}
		
		
		if(IntClust[i]=="CECc"){
			
			IntClustering[[i]]=CECc(List,distmeasure,normalize,method,t,r,nrclusters=nrclustersseq,weight,clust,linkage,alpha,Clustweight,StopRange)
			names(IntClustering)[i]="CECc"
			message("CECb")
		}
		}
		if(IntClust[i]=="WonM"){
			IntClustering[[i]]=WonM(List,type,distmeasure,normalize,method,nrclusters=nrclustersseq,clust,linkage,alpha,StopRange)				
			names(IntClustering)[i]="WonM"
			message("WonM")
		}
		
		if(IntClust[i]=="WeightedClust"){
			
			IntClustering[[i]]=WeightedClust(List,type,distmeasure,normalize,method,weight,Clustweight,clust,linkage,alpha,StopRange)
			names(IntClustering)[i]="WeightedClust"
			message("Weighted Clustering")
		}	
		
		if(IntClust[i]=="WeightedSim"){
			IntClustering[[i]]=WeightedSimClust(List,type,weight,clust=clust,linkage=linkage,alpha,distmeasure=distmeasure,normalize,method,gap=FALSE,maxK=50,nrclusters=nrclusters,names=c("Data1","Data2"),AllClusters=TRUE,StopRange)
			
			names(IntClustering)[i]="WeightedSim"
			message("WeightedSim")
		}
	}
	
	if(PlotCompare==TRUE){
		
		#Weights for CECa
		if(!(is.null(IntClustering$"CECa"))){
			grDevices::dev.new()
			ComparePlot(list(IntClustering$"CECa"),nrclusters,cols,fusionsLog,names=c(seq(1,0,-0.1)),main="CECa: Weights",...)
		}
		#Weights for CECb
		if(!(is.null(IntClustering$"CECb"))){
			grDevices::dev.new()
			ComparePlot(list(IntClustering$"CECb"),nrclusters,cols,fusionsLog,names=c(seq(1,0,-0.1)),main="CECb: Weights",...)
		}
		#Weights for CECc
		if(!(is.null(IntClustering$"CECc"))){
			grDevices::dev.new()
			ComparePlot(list(IntClustering$"CECc"),nrclusters,cols,fusionsLog,names=c(seq(1,0,-0.1)),main="CECc: Weights",...)
		}
		#Weights for WeightedClust
		if(!(is.null(IntClustering$"WeightedClust"))){
			grDevices::dev.new()
			ComparePlot(list(IntClustering$"WeightedClust"),nrclusters,cols,fusionsLog,names=c(seq(1,0,-0.1)),main="WeightedClust: Weights",...)
		}
		#Weights for WeightedSim
		if(!(is.null(IntClustering$"WeightedSim"))){
			grDevices::dev.new()
			ComparePlot(list(IntClustering$"WeightedSim"),nrclusters,cols,fusionsLog,names=c(seq(1,0,-0.1)),main="WeightedSim: Weights",...)
		}
		
		#All
		grDevices::dev.new()
		ComparePlot(IntClustering,nrclusters,cols,fusionsLog,WeightClust=TRUE,names=names(IntClustering),main="All Methods",...)
		
		grDevices::dev.new()
		L=list()
		L[[1]]=SeparateClustering[[1]]
		names(L)[[1]]="Data 1"
		for(k in 1:length(IntClustering)){
			L[[k+1]]=IntClustering[[k]]
			names(L)[k+1]=names(IntClustering)[k]
		}
		L[[length(L)+1]]=SeparateClustering[[2]]
		names(L)[[length(L)]]="Data 2"
		ComparePlot(L,nrclusters,cols,fusionsLog,WeightClust=TRUE,names=names(L),main="All Methods",...)
		
	}	
	out=c(SeparateClustering[1],SeparateClustering[2],IntClustering)
	return(out)
}
