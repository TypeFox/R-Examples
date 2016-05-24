clusterSumofSquares<- function(clusterOutput, sampledClusteringData){
   wss <- NULL
  for(k in 1:length(clusterOutput)){
   clusteredDatak <-  as.matrix(clusterOutput[[k]])
   wss[k] <- sum(apply(clusteredDatak,2,function(x) sum(rowSums(SumofSquares(arrayData=sampledClusteringData[x,])))))
  }
  return(wss)
}


SumofSquares<- function(arrayData){
  if(is.null(dim(arrayData))){ arrayData <- t(as.data.frame(arrayData))}
  overallMean <- sum(rowSums(arrayData))  /( nrow(arrayData) * ncol(arrayData) )
  overallMeanMat <- matrix(overallMean ,nrow = nrow(arrayData),ncol = ncol(arrayData))

  geneMeans <- rowSums(arrayData) / ncol(arrayData)
  geneMeansMat <- matrix(geneMeans ,nrow = nrow(arrayData),ncol = ncol(arrayData))

  doseMeans <- colSums(arrayData) / nrow(arrayData)
  doseMeansMat <- t(matrix(t(doseMeans ),ncol = nrow(arrayData),nrow = ncol(arrayData)))
  geneDoseResiduals <- (arrayData - geneMeansMat  - doseMeansMat
      + overallMeanMat)
  return(geneDoseResiduals^2)
}


clusterSumofAbsValues <- function(clusterOutput, sampledClusteringData){
   wss <- NULL
  for(k in 1:length(clusterOutput)){
   clusteredDatak <-  as.matrix(clusterOutput[[k]])
   wss[k] <- sum(apply(clusteredDatak,2,function(x) sum(rowSums(SumofAbsValues(arrayData=sampledClusteringData[x,])))))
  }
  return(wss)
}


SumofAbsValues <- function(arrayData){
geneDoseResiduals <- abs(medpolish(arrayData,trace.iter=FALSE)$residuals)
return(abs(geneDoseResiduals^2))
}


resampleORCME <- function(clusteringData, lambdaVector, robust=FALSE){
  nGene=100
  nData=100
  iter <- 1
  output <- list()
  while(iter <= nData){
    set.seed(12346*iter)
    ngene <- c(1:dim(clusteringData)[1])	
    sampledGenes <- sample(ngene,nGene,replace=TRUE)
    sampledClusteringData <- clusteringData[sampledGenes,]
    clusterOutput <- lapply(lambdaVector,function(x)ORCME(DRdata=sampledClusteringData,lambda=x,phi=2,robust=robust))
    nc <- sapply(c(1:length(clusterOutput)),function(x)dim(as.matrix(clusterOutput[[x]]))[2])
    
    if (robust==FALSE){    
	TSS <- sum(rowSums(SumofSquares(arrayData=sampledClusteringData)))
    	WSS <- clusterSumofSquares(clusterOutput=clusterOutput, sampledClusteringData=sampledClusteringData)
	} else {
	TSS <- sum(rowSums(SumofAbsValues(arrayData=sampledClusteringData)))
	WSS <- clusterSumofAbsValues(clusterOutput=clusterOutput, sampledClusteringData=sampledClusteringData)
	}

    output[[iter]] <-  cbind(lambda=lambdaVector,WSS=WSS,TSS=TSS,nc=nc)
    iter <- iter+1
  }
  return(output)
}

