SNFc<-function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,NN=20,mu=0.5,T=20,clust="agnes",linkage="ward",alpha=0.625,StopRange=FALSE){

	#Checking required data types and methods:
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	if(mu<0.3 | mu >0.8){
		message("Warning: mu is recommended to be between 0.3 and 0.8 for the SNF method. Default is 0.5.")
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
	
	
	#STEP 1: Distance Matrices
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		DistM=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize,method))
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		DistM=List
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
	}
	else{
		DistM=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
		
		OrderNames=rownames(DistM[[1]])
		for(i in 1:length(DistM)){
			DistM[[i]]=DistM[[i]][OrderNames,OrderNames]
		}
	}
	
	
	#STEP 2: Affinity Matrices
	
	AffinityMatrix.2=function (Diff, K = 20, sigma = 0.5) 
	{
		N = nrow(Diff)
		#Diff = (Diff + t(Diff))/2 #Why is this here?
		diag(Diff) = 0
		sortedColumns = as.matrix(t(apply(Diff, 2, sort)))
		finiteMean <- function(x) {
			mean(x[is.finite(x)])
		}
		means = apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
				.Machine$double.eps
		sum1 <- function(x, y) (x + y)
		Sig = outer(means, means, sum1)/3  + Diff/3 + .Machine$double.eps #Why times 2 and average if affinityMatrix?
		Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
		densities = stats::dnorm(Diff, 0, sigma * Sig, log = FALSE)
		densities=densities*(0.5*Sig*sqrt(2*pi))  #Rescale them back to have 1 on the diagonal? Or not necessary?
		W = (densities + t(densities))/2 #Why?
		W=densities
		return(W)	
	}
	
	AffM=lapply(seq(length(List)), function(i) AffinityMatrix.2(DistM[[i]], NN, mu))
	
	
	#STEP 3: Fuse Networks Into 1 Single Network
	#P and S matrix
	normalize <- function(X) {
		NMatrix=matrix(0,dim(X)[1],dim(X)[2])
		for(i in 1:dim(X)[1]){
			row=X[i,]
			row[i]=0
			D=sum(row)
			for(j in 1:dim(X)[2]){
				N=X[i,j]
			
				NMatrix[i,j]=N/(2*D)
			
				if(i==j){
					NMatrix[i,j]=1/2
				}
			}	
		
		}
		return(NMatrix)
	}	

	PMatrix= lapply(seq(length(List)),function(i) normalize(AffM[[i]]))	
	Wall=lapply(seq(length(List)),function(i) (AffM[[i]]+t(AffM[[i]]))/2)

	zero <- function(x,NN) {  #After affinityMatrix: the closest obs have the highest values.
		s = sort(x, index.return = TRUE)
		x[s$ix[1:(length(x) - NN)]] = 0
		return(x)
	}

	SMatrix=lapply(seq(length(List)),function(i) apply(AffM[[i]],1,zero, NN=NN))
	SMatrix=lapply(seq(length(List)),function(i) normalize(SMatrix[[i]])*2)

	nextW <- vector("list", length(List))
	for (i in 1:T) {	
	for (j in 1:length(List)) {
			sumWJ = matrix(0, dim(PMatrix[[j]])[1], dim(PMatrix[[j]])[2])
			for (k in 1:length(List)) {
				if (k != j) {
					sumWJ = sumWJ + PMatrix[[k]]
				}
			}
			nextW[[j]] = SMatrix[[j]] %*% (sumWJ/(length(List) - 1)) %*% t(SMatrix[[j]]) #update PMatrix
		}
		for (j in 1:length(List)) {
			PMatrix[[j]] = nextW[[j]] #+ diag(nrow(Wall[[j]]))  #why +1 on diagonals?
			Wall[[j]] = (Wall[[j]] + t(Wall[[j]]))/2           #why is this necessary?
			PMatrix[[j]]=normalize(PMatrix[[j]])   #normalize after every iteration
		}
	}

	SNF_FusedM=Reduce("+",PMatrix)
	SNF_FusedM=SNF_FusedM/length(List)
	SNF_FusedM = normalize(SNF_FusedM) #again normalization?
	SNF_FusedM = (SNF_FusedM + t(SNF_FusedM) + diag(nrow(SNF_FusedM)))/2

	rownames(SNF_FusedM)=rownames(List[[1]])
	colnames(SNF_FusedM)=rownames(List[[1]])


	Dist=1-SNF_FusedM
	
	#STEP 4: Perform Hierarchical Clustering with WARD Link
	
	HClust = cluster::agnes(Dist,diss=TRUE,method=linkage,par.method=alpha)
		
	
	#Output= list with the fused matrix and the performed clustering
	out=list(SNF_FusedM=SNF_FusedM,DistM=Dist,Clust=HClust)
	attr(out,'method')<-'SNF'
	return(out)
}

