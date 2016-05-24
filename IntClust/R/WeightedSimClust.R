WeightedSimClust<-function(List,type=c("data","dist","clusters"),weight=seq(0,1,0.01),clust="agnes",linkage=c("ward","flexible"),alpha=0.625,distmeasure=c("euclidean","tanimoto"),normalize=FALSE,method=NULL,gap=FALSE,maxK=50,nrclusters=NULL,names=c("B","FP"),AllClusters=FALSE,StopRange=FALSE,plottype="new",location=NULL){
	
	type<-match.arg(type)
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	if(class(weight)=="numeric" & length(weight)==length(List)){
		if(type=="data"){
			
			OrderNames=rownames(List[[1]])
			for(i in 1:length(List)){
				List[[i]]=List[[i]][OrderNames,]
			}
			#If given data matrices.
			##1: clustering on separate sources.
			##2: extract distance matrices.
			##3: determine nrclusters if not given.
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
			
			
			Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type=type,distmeasure[i],normalize=normalize,method=method,clust=clust,linkage=linkage[i],alpha=alpha,gap=gap,maxK=maxK,StopRange=StopRange))		
			
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
	
		}
		else if(type=="dist"){
			
			OrderNames=rownames(List[[1]])
			for(i in 1:length(List)){
				List[[i]]=List[[i]][OrderNames,OrderNames]
			}
			
			Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize,method,clust,linkage[i],alpha,gap,maxK,StopRange))
			
			Dist=List
			Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))		
			message("Warning: Are the distace matrices of the same range? If not, standardization is recommended.")	
		
		}
		else{
			Clusterings=List
			Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
			Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
			
			OrderNames=rownames(Dist[[1]])
			for(i in 1:length(Dist)){
				Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
			}
		}
		Weight=weight
		
	}
	else{
	
		temp=DetermineWeight_SimClust(List,type,weight,nrclusters,distmeasure,normalize,method,clust,linkage,alpha,gap,maxK,names,StopRange,plottype,location)

		Dist=lapply(seq(length(List)),function(i) temp[[1]][[i]]$DistM)
		
		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
			
		Weight=temp$Weight		
	}
	
	#Weighted Distance matrix
	
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)	
		return(temp)
	}
	DistW=weightedcomb(Weight,Dist)
	
	#Clustering on weighted distance matrix
		
	WeightedSimCluster=cluster::agnes(DistW,diss=TRUE,method="ward")

	out=list(Dist,Weight,DistW,Clust=WeightedSimCluster)	
	names(out)=c("Single_Distances","Weight","DistM","Clust")
	
	AllCluster=list()
	if(AllClusters==TRUE){
		if(class(weight)!="list"){
			condition<-function(l){		
				l=as.numeric(l)
				if( sum(l)==1 ){  #working with characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
					#return(row.match(l,t1))
					return(l)
				}
				else(return(0))
			}
			
			if(all(seq(1,0,-0.1)!=weight)){
				for(i in 1:length(weight)){
					rest=1-weight[i]
					if(!(rest%in%weight)){
						weight=c(weight,rest)
					}
				}
			}
			
			t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
			t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) 
			t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
			t4=t2[which(t3!=0)]
			weight=lapply(seq(length(t4)),function(i) rev(t4[[i]]))
			
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
			return(temp)
		}
		DistM=lapply(weight,weightedcomb,Dist)
		
		hclustW=lapply(seq(length(weight)),function(i) cluster::agnes(DistM[[i]],diss=TRUE,method="ward"))
		
		namesweights=c()
		for(i in 1:length(weight)){
			namesweights=c(namesweights,paste("Weight",weight[i],sep=" "))
		}
		Results=lapply(seq(1,length(hclustW)),function(i) return(c("DistM"=DistM[i],"Clust"=hclustW[i])))
		names(Results)=namesweights
				
		out=list(Dist, Weight,Results=Results,Clust=c("DistM"=DistW,"Clust"=WeightedSimCluster))	
		names(out)=c("Single Distances","Weight","Results","Clust")
	}	
	attr(out,'method')<-'WeightedSim'
	return(out)
}
