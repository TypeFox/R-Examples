
######################################################
# Soft Thresholding Funtion
######################################################
Soft <- function(a,b){
  
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a)*pmax(0,abs(a)-b))

}

######################################################
# Update the mean for different biclusters
######################################################
UpdateMus <- function(x, Cs, Ds, lambda=0){

  uniqCs <- sort(unique(Cs))
  uniqDs <- sort(unique(Ds))
  mus <- matrix(NA, nrow=length(uniqCs), ncol=length(uniqDs))

  for(k in uniqCs){
    for(r in uniqDs){
      if(lambda==0)  mus[k,r] <- mean(x[Cs==k,Ds==r])
      if(lambda>0) mus[k,r] <- Soft(mean(x[Cs==k,Ds==r]), lambda/(sum(Cs==k)*sum(Ds==r)))
      if(lambda<0) stop("Cannot have a negative tuning parameter value.")
    }
  }
  
  return(mus)

}

######################################################
# Update the cluster assignments
######################################################
UpdateClusters <- function(x, mus, curCs,curDs){
	
  Cs.new <- rep(NA, length(curCs))
  uniq <- 1:max(curCs)
  mus.expandcolumns<-mus[,curDs,drop=FALSE]

  for(i in 1:nrow(x)){
    dist2.clust <- NULL
    for(k in 1:length(uniq)){
      dist2.clust <- c(dist2.clust, sum((x[i,,drop=FALSE]-mus.expandcolumns[k,,drop=FALSE])^2))
    }
    wh <- which(dist2.clust==min(dist2.clust))
    Cs.new[i] <- wh[1]
  }
  
  return(Cs.new)
  
}


######################################################
# Renumber the cluster assignments to 1-4 in case cluster 
# number one is merged with some other clusters.
######################################################
ReNumber <- function(Cs){

  newCs <- rep(NA, length(Cs))
  uniq <- sort(unique(Cs))
  for(i in 1:length(uniq)){
    newCs[Cs==uniq[i]] <- i
  }
  return(newCs)

}

######################################################
# Objective function that we want to minimize
# Equation (2) in the paper
######################################################
Objective <- function(x, mus, Cs, Ds,lambda=0){

  return(sum((x-mus[Cs,Ds])^2)+2*lambda*sum(abs(mus)))

}


############################################
# Function to calculate BIC as in Section 5.2
############################################
CalculateBIC <- function(x,biclustobj){
	mat <- matrix(0, nrow=nrow(x), ncol=ncol(x))
	Cs <- biclustobj$Cs
	Ds <- biclustobj$Ds
	for(i in unique(Cs)){
		for(j in unique(Ds)){
			if(biclustobj$Mus[i,j]!=0){
				 mat[Cs==i,Ds==j] <- mean(x[Cs==i,Ds==j])
			}
		}
	}
	mat[biclustobj$mus==0] <- mean(x[biclustobj$mus==0])
	return(log(sum((x-mat)^2))*nrow(x)*ncol(x) + log(nrow(x)*ncol(x))*sum(biclustobj$Mus!=0))
}

#################################################
## Update Sigma using the graphical lasso
#################################################
updateSigma<-function(Delta,mus,x,alpha,Cs,Ds){
  Matrix<-matrix(0,nrow=length(Cs),ncol=length(Cs)) 
  uniqueCs<-sort(unique(Cs))
  uniqueDs<-sort(unique(Ds))
  
  for(i in uniqueCs){
  	S<-matrix(0,nrow=sum(Cs==i),ncol=sum(Cs==i)) 
  	for(j in uniqueDs){
      temp <- (matrix(x[Cs==i,Ds==j]-mus[Cs==i,Ds==j],nrow=sum(Cs==i),ncol=sum(Ds==j)))%*%solve(Delta[Ds==j,Ds==j])%*%t(matrix(x[Cs==i,Ds==j]-mus[Cs==i,Ds==j],nrow=sum(Cs==i),ncol=sum(Ds==j)))
      S <- S+temp
    } 
    S<-S/length(Ds)
    tempCovariance<-glasso(S,alpha,penalize.diagonal=TRUE)$w
    Matrix[Cs==i,Cs==i]<-tempCovariance
  }
  return(Matrix)
}


################################
## Update Delta using the graphical lasso
################################
updateDelta<-function(Sigma,mus,x,beta,Cs,Ds){
  Matrix<-matrix(0,nrow=length(Ds),ncol=length(Ds)) 
  uniqueCs<-sort(unique(Cs))
  uniqueDs<-sort(unique(Ds))
  
  for(i in uniqueDs){
#  	S<-matrix(0,nrow=ncol(x[,Ds==i]),ncol=ncol(x[,Ds==i])) 
  	S<-matrix(0,nrow=sum(Ds==i),ncol=sum(Ds==i)) 
  	for(j in uniqueCs){
      temp <- t(matrix(x[Cs==j,Ds==i]-mus[Cs==j,Ds==i],nrow=sum(Cs==j),ncol=sum(Ds==i)))%*%solve(Sigma[Cs==j,Cs==j])%*%(matrix(x[Cs==j,Ds==i]-mus[Cs==j,Ds==i],nrow=sum(Cs==j),ncol=sum(Ds==i)))
      S <- S+temp
    } 
    S<-S/length(Cs)
    tempCovariance<-glasso(S,beta,penalize.diagonal=TRUE)$w
    Matrix[Ds==i,Ds==i]<-tempCovariance
  }
  return(Matrix)
}


#############################################
## Update the mean matrix of MVN-biclustering
#############################################
updateMusMatrix<-function(x,Cs,Ds,lambda=0,Sigma,Delta){
  
  uniqueCs<-sort(unique(Cs))
  uniqueDs<-sort(unique(Ds))
  mus<-matrix(NA,nrow=length(uniqueCs),ncol=length(uniqueDs))	
  for(i in uniqueCs){
  	for(j in uniqueDs){
  	  invSigmak<-solve(Sigma[Cs==i,Cs==i])
  	  invDeltar<-solve(Delta[Ds==j,Ds==j])
  	  tempones<-matrix(1,nrow=sum(Cs==i),ncol=sum(Ds==j))
  	  denominator<-sum(diag(invSigmak%*%tempones%*%invDeltar%*%t(tempones)))
  	  numerator<-sum(diag(invSigmak%*%tempones%*%invDeltar%*%t(matrix(x[Cs==i,Ds==j],nrow=sum(Cs==i),ncol=sum(Ds==j)))))
  	  
  	  if(lambda==0) mus[i,j]<-numerator/denominator
  	  if(lambda>0) mus[i,j]<-Soft(numerator/denominator,lambda/denominator)
  	  if(lambda<0) stop("Negative tuning parameter is bad!")
  	}
  }

  return(mus)
}

################################
## Update Row Clusters
################################
UpdateRowCluster<-function(x,Sigma,Delta,mus,Cs,Ds){
  Cs.new<-rep(NA,length(Cs))	
  uniqueCs<-1:max(Cs)
  uniquemus<-unique(mus)
  invDelta<-solve(Delta)
  determinant<-sum(log(eigen(Delta,symmetric=TRUE,only.values=TRUE)$values))
  
  for(i in 1:nrow(x)){
  	dist2.clust<-NULL
  	
    for(j in 1:nrow(uniquemus)){
      temp<-(-0.5)*determinant-0.5*nrow(Delta)*log(Sigma[i,i])-0.5*t(x[i,]-uniquemus[j,])%*%invDelta%*%(x[i,]-uniquemus[j,])/Sigma[i,i]    
  	  dist2.clust<-c(dist2.clust,temp)
  	}
  	
  	wh<-which(dist2.clust==max(dist2.clust))
  	Cs.new[i]<-wh[1]
  }  	
  return(Cs.new)
}


################################
## Update Column Clusters
################################
UpdateColumnCluster<-function(x,Sigma,Delta,mus,Cs,Ds){
  Ds.new<-rep(NA,length(Ds))	
  uniqueDs<-1:max(Ds)
  uniquemus<-unique(t(mus))
  invSigma<-solve(Sigma)
  determinant<-sum(log(eigen(Sigma,symmetric=TRUE,only.values=TRUE)$values))
  
  for(i in 1:ncol(x)){
  	dist2.clust<-NULL
  	
  	for(j in 1:nrow(uniquemus)){
      temp<-(-0.5)*determinant-0.5*nrow(Sigma)*log(Delta[i,i])-0.5*t(x[,i]-uniquemus[j,])%*%invSigma%*%(x[,i]-uniquemus[j,])/Delta[i,i]
  	  dist2.clust<-c(dist2.clust,temp)
  	}
  	
  	wh<-which(dist2.clust==max(dist2.clust))
  	Ds.new[i]<-wh[1]
  }  	
  return(Ds.new)
}

###############################################
## Objective for calculating the log likelihood 
###############################################
MatrixObjective<-function(x,mus,Cs,Ds,Sigma,Delta,lambda=0,alpha=0,beta=0){
  tempSigma<-0
  tempSigma2<-0
  tempDelta<-0 
  tempDelta2<-0
  tempsum<-0 
  tempmus<-0
  
   
  for(i in 1:max(Cs)){
  	temp1<-log(det(solve(Sigma[Cs==i,Cs==i])))
  	tempSigma<-temp1+tempSigma
  	
  	temp11<-solve(Sigma[Cs==i,Cs==i])
  	tempSigma2<-sum(abs(temp11))+tempSigma2
  	}
  	
  for(j in 1:max(Ds)){
   	temp2<-log(det(solve(Delta[Ds==j,Ds==j])))
  	tempDelta<-temp2+tempDelta
  	
    temp21<-solve(Delta[Ds==j,Ds==j])
  	tempDelta2<-sum(abs(temp21))+tempDelta2
  	}	
  
  for(i in 1:max(Cs)){
  	for(j in 1:max(Ds)){
  		temp3<-solve(Sigma[Cs==i,Cs==i])%*%(matrix(x[Cs==i,Ds==j],nrow=sum(Cs==i),ncol=sum(Ds==j))-mus[Cs==i,Ds==j])%*%solve(Delta[Ds==j,Ds==j])%*%t(matrix(x[Cs==i,Ds==j],nrow=sum(Cs==i),ncol=sum(Ds==j))-mus[Cs==i,Ds==j])
  		tempsum<-sum(diag(temp3))+tempsum
  		
  		temp4<-abs(mus[Cs==i,Ds==j][1])
  		tempmus<-temp4+tempmus
  		}
  	}
  obj<--(ncol(x)/2)*tempSigma-(nrow(x)/2)*tempDelta+0.5*tempsum+lambda*tempmus+alpha*(ncol(x)/2)*tempSigma2+beta*(nrow(x)/2)*tempDelta2
  
  return(obj)
}

#########################################
## Function to deal with missing clusters
#########################################
ReNumberMatrix <- function(Cs){
  newCs <- rep(NA, length(Cs))
  uniq <- sort(unique(Cs))
  for(i in 1:length(uniq)){
    newCs[Cs==uniq[i]] <- i
  }
  return(newCs)
}

#####################################################
# Function to Calculate BIC for Matrix MVN
#####################################################
CalculateBICMatrix <- function(x,mclustering){
	mat <- matrix(0, nrow=nrow(x), ncol=ncol(x))
	Cs <- mclustering$Cs
	Ds <- mclustering$Ds
	for(i in unique(Cs)){
		for(j in unique(Ds)){
			if(mclustering$Mus[i,j]!=0){
				 mat[Cs==i,Ds==j] <- mean(x[Cs==i,Ds==j])
			}
		}
	}
	mat[mclustering$mus==0] <- mean(x[mclustering$mus==0])
	return(log(sum((x-mat)^2))*nrow(x)*ncol(x) + log(nrow(x)*ncol(x))*sum(mclustering$Mus!=0))
}

