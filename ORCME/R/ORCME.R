`addRows` <- function(originalData, clusteredData, robust) {

  HCluster <- Hcriterion(clusteredData, robust=robust)$`overallHvalue `

  if(robust==FALSE){
	  overallMeanCluster <- sum(clusteredData) / (nrow(clusteredData) * ncol(clusteredData))
	  doseMeansCluster <- colSums(clusteredData) / nrow(clusteredData)
	  newGenesMeans <- rowSums(originalData) / ncol(originalData)
	  list <- seq(1:nrow(originalData))
	  ResidualsData <- unlist(lapply(list, function (x) mean(((originalData[x,] - newGenesMeans[x]) - (doseMeansCluster - overallMeanCluster))^2)))

  } 
  if(robust==TRUE){
	  mp <- medpolish(clusteredData,trace.iter=FALSE)
	  overallMeanCluster <- mp$overall
	  doseMeansCluster <- mp$col
	  newGenesMeans <- medpolish(originalData,trace.iter=FALSE)$row
	  list <- seq(1:nrow(originalData))
	  ResidualsData <- unlist(lapply(list, function (x) mean(as.numeric(abs((originalData[x,] - newGenesMeans[x]) - (doseMeansCluster - overallMeanCluster))))))
  }

  geneAddition <- ResidualsData < HCluster
  geneAddition <- geneAddition[!(rownames(originalData)[geneAddition] %in% row.names(clusteredData))] 

    if(sum(geneAddition)!=0){
       addData <- originalData[geneAddition, ]
       #clusteredData$ID <- rownames(clusteredData)
       #addData$ID <- rownames(addData)
       clusteredOutput <-  rbind(clusteredData,addData)
       #clusteredOutput <-  clusteredOutput[!duplicated(clusteredOutput$ID),]
       #clusteredOutput <- clusteredOutput[,-5]
       dataLogic <-   row.names(originalData) %in% row.names(clusteredOutput)
       dataOutput <- originalData[!dataLogic, ]
  }else{
       clusteredOutput <-  clusteredData
       dataLogic <-   row.names(originalData) %in% row.names(clusteredData)
       dataOutput <- originalData[!dataLogic,]
  }
  output <- list(clusteredOutput,dataOutput)
  names(output) <- c("clusteredOutput","dataOutput")
  return(output)
}

'clusterLogic' <- function(incClusterData,pnames,x){
     idlogic <- pnames %in%  incClusterData[[x]]
     return(idlogic)
}

`clusteredMem` <- function(parentData, incClusterData){
    parentData <- parentData
    ncluster <- c(1:length(incClusterData))
    pnames <- row.names(parentData)
    clustoutput <- sapply(ncluster,function(x) clusterLogic(incClusterData,pnames,x))
    return(clustoutput)
}


`Hcriterion` <- function(arrayData, addRow=FALSE, robust){

	if(robust==FALSE){
  		overallMean <- sum(arrayData) / ( nrow(arrayData) * ncol(arrayData) )
  		overallMeanMat <- matrix(overallMean , nrow = nrow(arrayData), ncol = ncol(arrayData))

  		geneMeans <- rowSums(arrayData) / ncol(arrayData)
  		geneMeansMat <- matrix(geneMeans, nrow = nrow(arrayData),ncol = ncol(arrayData))

  		doseMeans <- colSums(arrayData) / nrow(arrayData)
  		doseMeansMat <- t(matrix(t(doseMeans ), ncol = nrow(arrayData),nrow = ncol(arrayData)))
  		geneDoseResiduals <- (arrayData - geneMeansMat  - doseMeansMat  + overallMeanMat)^2
  		geneResiduals <- rowSums(geneDoseResiduals) / ncol(geneDoseResiduals)
  		overallHvalue <- sum(geneDoseResiduals) / (nrow(geneDoseResiduals) * ncol(geneDoseResiduals))
	} 
	if(robust==TRUE){
		geneDoseResiduals <- abs(medpolish(arrayData,trace.iter=FALSE)$residuals)
		geneResiduals <- rowSums(geneDoseResiduals) / ncol(geneDoseResiduals)
  		overallHvalue <- sum(geneDoseResiduals) / (nrow(geneDoseResiduals) * ncol(geneDoseResiduals))
	}
	
  mseData <- list(geneResiduals ,overallHvalue )
  names(mseData) <- c("geneResiduals ","overallHvalue ")
  return(mseData)
}


`deltaInit` <- function(arrayData,lambda, robust){
 
  initialHvalue <- Hcriterion(arrayData, robust=robust)
  delta <- lambda * initialHvalue$overallHvalue  

  return(delta)
}

`singleNodeDeletion` <- function(arrayData,delta,phi, robust){
  sndData <-Hcriterion(arrayData, robust=robust)
  overallHvalue <- sndData$overallHvalue
  geneResiduals <- sndData$geneResiduals
  fuse<-1
    # the fuse: in case of the resampling with replacement, it could happen 
    # that the data set have following property:
    # In one point, the data set could reduce from few rows lines to zero rows
    # in one step and then the algorithm would collapse since the fuctions 
    # would not be defined => the fuse prevents this situation 

  while(overallHvalue > delta & dim(arrayData)[1]>phi & fuse>0)  {
      geneResidualsLogic <- geneResiduals < max(geneResiduals)
	  
      if(dim(arrayData[geneResidualsLogic , ])[1]==0){fuse <- 0
	} else{
	   arrayData <- arrayData[geneResidualsLogic , ]
       sndData <- Hcriterion(arrayData, robust=robust)
       overallHvalue <- sndData$overallHvalue
       geneResiduals <- sndData$geneResiduals
	 }
 }

  return(arrayData)

}



`ORCME` <- function(DRdata,lambda,phi, robust=FALSE){
  clusterRowNames <- list()
  DRdata <- as.data.frame(DRdata)
  row.names(DRdata) <- paste("g",c(1:dim(DRdata)[1]),sep="")
  DDdata <- DRdata
  delta <- deltaInit(arrayData = DDdata, lambda, robust=robust)
  phi <- phi
  iter <- 1
  while (dim(DDdata)[1] > phi){

    ## apply single node deletion algorithm of Cheng and Church, 2000
    clsData <- singleNodeDeletion(arrayData=DDdata, delta, phi, robust=robust)

    ## apply algorithm for row addition by Cheng and Church, 2000
    clusterAddRows <- addRows(originalData = DDdata, clusteredData =clsData, robust=robust)

    ## discovered cluster
    clusterRowNames[[iter]] <- row.names(clusterAddRows$clusteredOutput)
    ## remaining data to be clustered
    DDdata <- clusterAddRows$dataOutput
    iter <- iter+1
  }
  if(!is.null(dim(DDdata))){clusterRowNames[[iter]] <- row.names(DDdata) }
  clusout <- clusteredMem(parentData=DRdata,incClusterData=clusterRowNames)
  output <- as.data.frame(clusout[,colSums(clusout)!=0])
  row.names(output) <- row.names(DRdata)
  return(output)
}


EstimatedParameters <- function(DataSet, ORCMEoutput, clusterID, robust=FALSE){
	useData <- DataSet[ORCMEoutput[,clusterID],]
	if(robust==FALSE){
  		overallMean <- sum(useData) / ( nrow(useData) * ncol(useData) )
  		rowEffects <- rowSums(useData) / ncol(useData) - overallMean
  		colEffects <- colSums(useData) / nrow(useData) - overallMean
	} else {
		mp <- medpolish(useData,trace.iter=FALSE)
  		overallMean <- mp$overall
  		rowEffects <- mp$row
  		colEffects <- mp$col
	}
	result <- list(overallMean, rowEffects, colEffects)
	names(result) <- c("overallMean", "rowEffects", "colEffects")
return(result)
}
