WeightedClust <- function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,weight=seq(1,0,-0.1),WeightClust=0.5,clust="agnes",linkage="ward",alpha=0.625,StopRange=FALSE){ # weight = weight to data1

	
	#Step 1: compute distance matrices:
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
	
	
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize,method))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
	}
	else{
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$Dist))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		OrderNames=rownames(DistM[[1]])
		for(i in 1:length(DistM)){
			DistM[[i]]=DistM[[i]][OrderNames,OrderNames]
		}
	}
	
	#Step 2: Weighted linear combination of the distance matrices:
	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))
	
	}
	else if(class(weight)=='list' & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	
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
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		t4=t2[which(t3!=0)]
		weight=lapply(seq(length(t4)),function(i) rev(t4[[i]]))
		
	}
	if(class(weight)=="list" & "x" %in% weight[[1]]){ #x indicates a free weight
		newweight=list()
		for(k in 1:length(weight)){
			w=weight[[k]]
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
			
			newweight[k]=new
		}
		
		weight=newweight
	}
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)	
		return(temp)
	}
	
	DistClust=NULL
	Clust=NULL
	
	DistM=lapply(seq(length(weight)),function(i) weightedcomb(weight[[i]],Dist=Dist))
	namesweights=c()
	WeightedClust=lapply(seq(length(weight)),function(i) cluster::agnes(DistM[[i]],diss=TRUE,method=linkage,par.method=alpha))
	for(i in 1:length(WeightedClust)){
		namesweights=c(namesweights,paste("Weight",weight[i],sep=" "))
		if(all(weight[[i]]==WeightClust)){
			Clust=WeightedClust[[i]]	
			DistClust=DistM[[i]]
		}
	}	
	
	if(is.null(DistClust)){
		DistClust=weightedcomb(WeightClust,Dist=Dist)
		Temp=cluster::agnes(DistClust,diss=TRUE,method=linkage,par.method=alpha)
		Clust=Temp
	}
	
	Results=lapply(seq(1,length(WeightedClust)),function(i) return(c("DistM"=DistM[i],"Clust"=WeightedClust[i])))
	names(Results)=namesweights
	
	# return list with objects
	out=list(Dist=Dist,Results=Results,Clust=list("DistM"=DistClust,"Clust"=Clust))
	attr(out,'method')<-'Weighted'
	return(out)
	
}
