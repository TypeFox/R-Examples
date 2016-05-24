DetermineWeight_SimClust<-function(List,type=c("data","dist","clusters"),weight=seq(0,1,by=0.01),nrclusters=NULL,distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,clust="agnes",linkage=c("ward","ward"),alpha=0.625,gap=FALSE,maxK=50,names=c("B","FP"),StopRange=FALSE,plottype="new",location=NULL){
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	type<-match.arg(type)
	if(type=="data"){
		#for(a in 1:length(distmeasure)){
		#	if(distmeasure[[a]]=='euclidean'){
		#		stand<-function(c){
		#			minc=min(c)
		#			maxc=max(c)
		#			c1=(c-minc)/(maxc-minc)
		#			return(c1)				
		#		}
		#		List[[a]]=apply(List[[a]],2,stand)			
		#	}
		#	
		#}
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize,method,clust,linkage[i],alpha,gap,maxK,StopRange))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		out<-list(Clusterings)
		names(out)="ClusterSep"		
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize,method,clust,linkage[i],alpha,gap,maxK,StopRange))
		
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		out<-list(Clusterings)
		names(out)="ClusterSep"		
		
	}
	else{
		
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		out<-list(Clusterings)
		names(out)="ClusterSep"		
		
	}
	
	
	namesw=c()
	for(i in 1:length(names)){
		namesw=c(namesw,paste("w_",names[i],sep=""))
	}
	namesJ=c()
	for(i in 1:length(names)){
		namesJ=c(namesJ,paste("J(sim",names[i],",simW)",sep=""))
	}
	namesR=c()
	combs=utils::combn(seq(length(List)),m=2,simplify=FALSE)
	for(i in 1:length(combs)){
		namesR=c(namesR,paste("J_",names[combs[[i]][1]],"/J_",names[combs[[i]]][2],sep=""))
	}
	
	labels<-c(namesw,namesJ,namesR)
	
	ResultsWeight<-matrix(0,ncol=length(labels),nrow=length(weight))
	#data.frame(col1=numeric(),col2=numeric(),col3=numeric(),col4=numeric())
	colnames(ResultsWeight)=labels
	
	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))		
	}
	else if(class(weight)=='list' & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	else{
		message('The weights are considered to be a sequence, each situation is investigated')
	}
	
	if(class(weight)!="list"){
		condition<-function(l){		
			l=as.numeric(l)
			if( sum(l)==1 ){  #working wit characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		weight=t2[which(t3!=0)]
	}
	
	if(class(weight)=="list" & "x" %in% weight[[1]]){ #x indicates a free weight
		for(i in 1:length(weight)){
			w=weight[[i]]
			weightsfordata=which(w!="x") #position of the provided weight = position of the data to which the weight is given
			givenweights=as.numeric(w[weightsfordata])
			
			stilltodistribute=1-sum(givenweights)
			
			newweights=seq(stilltodistribute,0,-0.1)
			
			t1=gtools::permutations(n=length(newweights),r=length(List)-length(weightsfordata),v=as.character(newweights),repeats.allowed = TRUE)
			Input1=as.list(seq_len(nrow(t1)))
			Input2=lapply(seq(length(Input1)),function(i) {Input1[[i]][length(Input1[[i]])+1]=stilltodistribute
														   return(Input1[[i]])})
			t2=lapply(seq(length(Input2)), FUN=function(i){if(sum(as.numeric(t1[Input2[[i]][1],])+0.00000000000000002775)==Input2[[i]][2]) return(as.numeric(t1[i,])) else return(0)}) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
			t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
			weightsforotherdata=t2[which(t3!=0)]
			
			new=list()
			for(i in 1:length(weightsforotherdata)){
				w1=weightsforotherdata[[i]]
				new[[i]]=rep(0,length(List))
				new[[i]][weightsfordata]=givenweights
				new[[i]][which(new[[i]]==0)]=w1
			}
			
			weight=new
		}
	}
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)		
	}
	DistM=lapply(weight,weightedcomb,Dist)
	
	hclustOr=lapply(seq(length(List)),function(i) stats::cutree(Clusterings[[i]]$Clust,nrclusters))
	hclustW=lapply(seq(length(weight)),function(i) stats::cutree(cluster::agnes(DistM[[i]],diss=TRUE,method="ward"),nrclusters))
	
	Counts=function(clusterlabs1,clusterlabs2){
		index=c(1:length(clusterlabs1))
		allpairs=utils::combn(index,2,simplify=FALSE)  #all pairs of indices: now check clutserlabels for every pair==> only 1 for loop
		n11=n10=n01=n00=0
		
		counts<-function(pair){
			if(clusterlabs1[pair[1]]==clusterlabs1[pair[2]]){
				if(clusterlabs2[pair[1]]==clusterlabs2[pair[2]]){
					n11=n11+1
				}
				else{
					n10=n10+1						
				}
			}
			else{
				if(clusterlabs2[pair[1]]==clusterlabs2[pair[2]]){
					n01=n01+1
				}
				else{
					n00=n00+1
				}
				
			}
			return(c(n11,n10,n01,n00))		
		}
		
		n=lapply(seq(length(allpairs)),function(i) counts(allpairs[[i]]))
		nn=Reduce("+",n)
		#2: compute jaccard coefficient	
		Jac=nn[1]/(nn[1]+nn[2]+nn[3])
		return(Jac)
	}
	
	Jaccards<-function(hclust){
		jacs=lapply(seq(length(hclustOr)),function(i) Counts(clusterlabs1=hclustOr[[i]],clusterlabs2=hclust))
		return(unlist(jacs))
	}
	
	AllJacs=lapply(hclustW,Jaccards)  #make this faster:lapply + transfrom to data frame with plyr package
	
	Ratios<-function(Jacs){	
		combs=utils::combn(seq(length(List)),m=2,simplify=FALSE)
		ratio<-function(v,Jacs){
			return(Jacs[v[1]]/Jacs[v[2]])
		}
		
		ratios=lapply(seq(length(combs)),function(i) ratio(v=combs[[i]],Jacs=Jacs))
		
	}
	
	AllRatios=lapply(seq(length(AllJacs)),function(i) unlist(Ratios(AllJacs[[i]])))

	
	ResultsWeight=t(mapply(c,weight,AllJacs,AllRatios))
	colnames(ResultsWeight)=labels
	
	#Choose weight with ratio closest to one==> smallest where this happens: ##### START HERE WITH OPTIMIZATION #####
	ResultsWeight=cbind(ResultsWeight,rep(0,nrow(ResultsWeight)))
	colnames(ResultsWeight)[ncol(ResultsWeight)]="trick"
	Weight=ResultsWeight[which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),c(1:length(List))]
	
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
	graphics::plot(x=0,y=0,type="n",xlim=c(0,dim(ResultsWeight)[1]),ylim=c(min(ResultsWeight[,namesR]),max(ResultsWeight[,namesR])),xlab="",ylab="Ratios")
	if(is.null(ncol(ResultsWeight[,namesR]))){
		L=1
	}
	else{
		L=ncol(ResultsWeight[,namesR])
	}
	for(i in 1:L){
		graphics::lines(1:dim(ResultsWeight)[1],y=ResultsWeight[,namesR[i]],col=i)
	}
	graphics::abline(h=0,v=which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),col="black",lwd=2)
	graphics::mtext("Weight Combinations", side=1, line=3)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),line=1,tck=1,lwd=2)
	plottypeout(plottype)
	
	ResultsWeight=ResultsWeight[,-ncol(ResultsWeight)]
	out[[2]]=ResultsWeight
	out[[3]]=Weight
	names(out)=c("ClusterSep","Result","Weight")
	
	
	
	return(out)
	
}
