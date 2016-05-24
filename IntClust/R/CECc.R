CECc<-function(List,distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,t=10,r=NULL,nrclusters=NULL,weight=NULL,clust="agnes",linkage=c("flexible","flexible"),alpha=0.625,WeightClust=0.5,StopRange=FALSE){
	
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
		
	if(is.null(nrclusters)){
		stop("Give a number of clusters to cut the dendrogram into for each data modality.")
	}
	

	#Step 1: Take a random subset of features from A1 and A2
	#Notation facility:
	
	#Put all data in the same order
	OrderNames=rownames(List[[1]])
	for(i in 1:length(List)){
		List[[i]]=List[[i]][OrderNames,]
	}
	
	#Put up Incidence matrix for each data modality
	nc=c()
	Incidence=list()
	for (i in 1:length(List)){
		Incidence[[i]]=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
		rownames(Incidence[[i]])=rownames(List[[i]])
		colnames(Incidence[[i]])=rownames(List[[i]])
		nc=c(nc,ncol(List[[i]]))
	}
	evenn=function(x){if(x%%2!=0)x=x-1 else x}
	nc=lapply(nc,FUN=evenn)
	nc=unlist(nc)
	
	#Repeat for t iterations
	for(g in 1:t){
		message(g)	
		if(is.null(r)){
			r=c()
			for(i in 1:length(nc)){
				r=c(r,sample((nc[i]/2):(nc[i]-1),1))
			}
		}
		
		#take random sample:
		A_prime=list()
		for(i in 1:length(r)){
			A=List[[i]]
			temp=sample(ncol(A),r[i],replace=FALSE)
			A_prime[[i]]=A[,temp]
			
			Ok=FALSE
			while(Ok==FALSE){
				if(class(as.vector(as.matrix(A_prime[[i]])))=="character"){
					if(distmeasure=="MCA_coord"){
						dat=t(A_prime[[i]])
						dat=as.data.frame(dat)
						dat=dat[, sapply(dat, nlevels) > 1,drop=FALSE]					
						if(ncol(dat)<=1){
							temp=sample(ncol(A),r[i],replace=FALSE)
							A_prime[[i]]=A[,temp]
						}
						else{
							A_prime[[i]]=t(dat)
							Ok=TRUE
						}
					}
				}
				
				else if(any(rowSums(A_prime[[i]])==0)){
					temp=sample(ncol(A),r[i],replace=FALSE)
					A_prime[[i]]=A[,temp]
				}
				else{
					Ok=TRUE
				}				
			}			
			
		}
		
		
			
		#Step 2: apply hierarchical clustering on A1_prime and A2_prime + cut tree into nrclusters
		
		DistM=lapply(seq(length(A_prime)),function(i) Distance(A_prime[[i]],distmeasure=distmeasure[i],normalize,method))
		
		CheckDist<-function(Dist,StopRange){
			if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
				message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				Dist=Normalization(Dist,method="Range")
			}
			else{
				Dist=Dist
			}
		}
		
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
		
		HClust_A_prime=lapply(seq(length(DistM)),function(i) cluster::agnes(DistM[[i]],diss=TRUE,method=linkage[i],par.method=alpha))
				
		for(k in 1:length(nrclusters)){
			MembersofClust=list()
			for (i in 1:length(HClust_A_prime)){
				Temp=stats::cutree(HClust_A_prime[[i]],nrclusters[k])	
				MembersofClust[[i]]=matrix(1,dim(List[[i]])[1],dim(List[[i]])[1])
				
				for(l in 1:length(Temp)){
					label=Temp[l]
					sameclust=which(Temp==label)
					MembersofClust[[i]][l,sameclust]=0
				}
				Incidence[[i]]=Incidence[[i]]+MembersofClust[[i]]
			}	
		}
		
		
	}
	
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
	IncidenceComb=lapply(weight,weightedcomb,Incidence)
	
	
	namesweights=c()
	CEC=list()
	for (i in 1:length(IncidenceComb)){
		CEC[[i]]=cluster::agnes(IncidenceComb[[i]],diss=TRUE,method="ward")
		namesweights=c(namesweights,paste("Weight",weight[i],sep=" "))
		if(all(weight[[i]]==WeightClust)){
			Clust=CEC[i]
			DistClust=IncidenceComb[i]
		}
	}
	
	Results=lapply(seq(1,length(CEC)),function(i) return(c("DistM"=IncidenceComb[i],"Clust"=CEC[i])))
	names(Results)=namesweights
	
	out=list(Incidence=Incidence,Results=Results,Clust=c("DistM"=DistClust,"Clust"=Clust))
	attr(out,'method')<-'CEC'
	return(out)	
	
}
