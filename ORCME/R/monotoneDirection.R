## monotonicity
minmax <- function(tmpMatrix){
  tmp<-rep(NA,ncol(tmpMatrix))
  for(i in 1:ncol(tmpMatrix)){
   tmp[i] <- max(tmpMatrix[,i])
  }
  return(min(tmp))
}


isotonicmeans <- function(mean.y,invsigma,isoDir){
k <- length(mean.y)
mu.data <- rep(NA,k)

if(isoDir == "up"){
  for(i in 1:k){
    L <- c(1:i)
    U <- c(i:k)
    tmpmat <- matrix(NA,nrow = length(L),ncol=length(U))
    for(j in 1:length(L)) {
      for(l in 1:length(U)){

        tmpidentity <- rep(1,(U[l]-L[j]+1))
        tmpinvsigma <- invsigma[L[j]:U[l],L[j]:U[l]]
        tmpmu <- mean.y[L[j]:U[l]]
        tmpnumerator <-  t(tmpidentity) %*%  tmpinvsigma %*% tmpmu
        tmpdenominator <- t(tmpidentity) %*%  tmpinvsigma %*% tmpidentity
        tmpmat[j,l] <-   tmpnumerator/ tmpdenominator
      }

    }
        mu.data[i] <- minmax(tmpMatrix=tmpmat)[1]
  }
 }else if(isoDir == "dn"){
for(i in k:1){
    L <- c(k:i)
    U <- c(i:1)
    tmpmat <- matrix(NA,nrow = length(L),ncol=length(U))
    for(j in 1:length(L)) {
      for(l in 1:length(U)){

        tmpidentity <- rep(1,(L[j]-U[l]+1))
        tmpinvsigma <- invsigma[L[j]:U[l],L[j]:U[l]]
        tmpmu <- mean.y[L[j]:U[l]]
        tmpnumerator <-  t(tmpidentity) %*%  tmpinvsigma %*% tmpmu
        tmpdenominator <- t(tmpidentity) %*%  tmpinvsigma %*% tmpidentity
        tmpmat[j,l] <-   tmpnumerator/ tmpdenominator
      }

    }
        mu.data[i] <- minmax(tmpMatrix=tmpmat)[1]
  }

}
  return( mu.data)
}

isominmax <- function (compData, doseData) 
{
	compData2 <- compData[, order(doseData, decreasing = F)]
	doseData2 <- doseData[order(doseData, decreasing = F)]
	isotmp <- rep(NA, nrow(compData2))
	isomeans <- matrix(NA, nrow = nrow(compData2), ncol = length(unique(doseData2)))
	for (i in 1:nrow(compData2)) {
		ndose <- tapply(doseData, doseData, length)
		exampleminmax <- as.vector(compData2[i, ])
		mean.exampleminmax <- tapply(exampleminmax, doseData2, mean)
		upcmean.minmax <- pava(mean.exampleminmax, 
				decreasing=FALSE, long.out=FALSE, stepfun=FALSE)
		downcmean.minmax <- pava(mean.exampleminmax, 
				decreasing=TRUE, long.out=FALSE, stepfun=FALSE)
		
		likelihoodmu.dose <- rep(mean.exampleminmax[1], ndose[1])
		upmu.dose <- rep(upcmean.minmax[1], ndose[1])
		downmu.dose <- rep(downcmean.minmax[1], ndose[1])
		for (j in 2:length(mean.exampleminmax)) {
			likelihoodmu.dose <- c(likelihoodmu.dose, rep(mean.exampleminmax[j], 
							ndose[j]))
			upmu.dose <- c(upmu.dose, rep(upcmean.minmax[j], 
							ndose[j]))
			downmu.dose <- c(downmu.dose, rep(downcmean.minmax[j], 
							ndose[j]))
		}
		nullvar <- sum((exampleminmax - likelihoodmu.dose)^2)
		upvar <- sum((exampleminmax - upmu.dose)^2)
		dnvar <- sum((exampleminmax - downmu.dose)^2)
		uplike <- (nullvar - upvar)
		dnlike <- (nullvar - dnvar)
		isotmp[i] <- c("up", "dn")[c(uplike, dnlike) == max(c(uplike, 
								dnlike))]
		if (isotmp[i] == "up") {
			isomeans[i, ] <- upcmean.minmax
		}  else {
			isomeans[i, ] <- downcmean.minmax
		}
	}
	outtmp <- list(isotmp, isomeans)
	names(outtmp) <- c("Direction", "isomeans")
	return(outtmp)
}



`monotoneDirection` <- function(geneData, doseData){
  geneDataMat <- as.matrix(geneData)
  H1 <- isominmax(compData=geneDataMat,doseData=doseData)
  gdir <- as.vector(H1$Direction)
  arrayMean <- H1$isomeans
  incData <- arrayMean [gdir == "up",]
  decData <- arrayMean [gdir == "dn",]
  geneincData <- geneData [gdir == "up",]
  genedecData <- geneData [gdir == "dn",]
  dirData <- list(gdir,incData,decData,geneincData ,genedecData,arrayMean)
  names(dirData)<-c('direction','incData','decData','obsincData','obsdecData','arrayMean')
  
  return(dirData)
}

