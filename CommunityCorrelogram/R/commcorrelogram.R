# TODO: Add comment
# 
# Author: Tim Kuehlhorn and Malia Andrus
###############################################################################
setClass('commcorrelogram',representation(community.correlogram='data.frame'),
	prototype=prototype(community.correlogram=data.frame(lag.distance=numeric(),Num.Pairs=numeric(),
	statistic=numeric(),significance=numeric())))

setValidity('commcorrelogram',
	function(object){
		dim(object@community.correlogram)[2]==4 && 
		((abs(object@community.correlogram[,3])<=1 && object@community.correlogram[,4]<=1) |
		 (is.na(object@community.correlogram[,3]) && is.na(object@community.correlogram[,4])))
	})

setMethod("plot",signature(x='commcorrelogram',y='missing'),
	definition=plot.commcorrelogram)

setGeneric('mod',def=function(object,...) object)

setMethod('mod',signature=signature(object='commcorrelogram'),
	definition=mod.commcorrelogram)


azimuthTest <- function(direction, aloc, bloc) {
    pairVector = as.numeric(aloc - bloc)
    xyNorm <- norm(matrix(pairVector[1:2], 1), "f")
    aziTest <- 0
    dipTest <- 0
    if(xyNorm != 0) {
        azimuthVector <- cbind(cos(direction$azimuth), sin(direction$azimuth))
        azimuthDot <- sum( (pairVector[1:2] / xyNorm) * azimuthVector)       
        if(round(abs(azimuthDot),10) >= round(abs(cos(direction$azimuthTol)),10)) {
            pairBandwidth <- xyNorm - abs(azimuthDot*xyNorm)
           if(abs(pairBandwidth) <= direction$bandwidth) {
               aziTest <- 1
           } 
        } 
       } else if (xyNorm == 0) {
        aziTest <- 1
    		}  
    if (aziTest == 1 & norm(matrix(pairVector,1), "f") != 0) {
        dipVector <- cbind(round(cos(direction$dipAngle),10), round(sin(direction$dipAngle),10))
        declination <- sum( c(xyNorm,pairVector[3])/norm(matrix(pairVector,1), "f") * dipVector)
    		if (declination >= round(abs(cos(direction$dipTol)),10)) {
                vBandwidth <- norm(matrix(pairVector,1), "f") - declination*norm(matrix(pairVector,1), "f")
                if( abs(vBandwidth) <= direction$dipBandwidth) {
                    dipTest <- 1
                }
            } 
    }else if (norm(matrix(pairVector,1), "f") == 0) {
    	    	dipTest <- 1
    }
    return(as.numeric(aziTest==1 & dipTest==1))
}

locationEquality <- function(pairNames) {
    as.numeric(pairNames[1] == pairNames[2])
}

distanceLags <- function(lagNumber, lagSize, lagTol, distance) {
    lags <- rep(0,lagNumber+1)
    if (distance == 0) lags[1] <- 1
    for(i in c(2:(lagNumber+1))) {
        spacing <- (i-2) * lagSize
		if(distance >= spacing - lagTol && distance < spacing + lagTol && distance!=0) {
            lags[i] <- 1
        }
    }
    return(lags)
}


RValues <- function(pairFrame,lagNumber,metric,mantmeth='spearman') {
	pairFrame<-data.frame(pairFrame)
    lagRValues <- rep(0,lagNumber+1)
    NumPairs <- rep(0,lagNumber+1)
    lagAvgDist <- rep(0,lagNumber+1)
    if (metric=='anosim') {
    	pairFrame[,3] <- rank(pairFrame[,3]) 
    	withoutTotal <- sum(pairFrame[,3])
    	}
    if (metric=='mantel'){
    	pairFrame[,3] <- scale(pairFrame[,3])
    }
    for (i in c(1:(lagNumber+1))) {  
        lagSet <- pairFrame[(pairFrame[,(3+i)]==1 & pairFrame[,1] == 1),]
        numPairs <- nrow(lagSet)
        if(numPairs != 0) {
    		if(metric=='anosim'){
            	lagWithinSum <- sum(lagSet[,3])
            	lagWithoutSum <- withoutTotal - lagWithinSum
            	lagWithinAvg <- lagWithinSum / nrow(lagSet)
            
            	if(nrow(pairFrame) - nrow(lagSet) != 0) {
                	lagWithoutAvg <- lagWithoutSum / (nrow(pairFrame) - nrow(lagSet))
            	} else {
                	lagWithoutAvg <- 0
            	}
            	lagRValues[i] <- (lagWithoutAvg - lagWithinAvg) / (nrow(pairFrame)/2)
    		}
    		else if(metric=='mantel'){
    			lagRValues[i] <- suppressWarnings(cor.test(as.numeric(!(pairFrame[,(3+i)]==1 & pairFrame[,1] == 1)),pairFrame[,3]
    				,method=mantmeth))$estimate   		
    			}
            
        } else {
            lagRValues[i] <- NA
        }
    NumPairs[i] <- numPairs
    lagAvgDist[i] <- mean(lagSet[,2])
    }
    return(cbind(lagRValues,NumPairs,lagAvgDist))
}

lag.sort<-function(sampleData, sampleTime=NULL, 		
	sampleLocation=NULL,LocationNames=NULL,option = 1,lagNumber, 
	lagSize, lagTol, anisotropic = FALSE, azimuth, 		
	azimuthTol,bandwidth,dipAngle,dipTol,dipBandwidth,distmeth='bray'){
	   pairEquality <- rep(1,nrow(sampleData)*(nrow(sampleData)-1)*0.5)
	    if (option %in% c(1,3)){
    	distanceMatrix <- vegdist(sampleLocation, na.rm=TRUE, method='euclidean')
    	distanceVector <- as.vector(distanceMatrix)
    	simMatrix <- vegdist(sampleData, method=distmeth, na.rm=TRUE) 
    	if (anisotropic == TRUE){
    		pairlist<-combn(c(1:dim(sampleLocation)[1]),2)
    		direction<-list(azimuth=azimuth*(pi/180),azimuthTol=azimuthTol*(pi/180),bandwidth=bandwidth,dipAngle=dipAngle*(pi/180),dipTol=dipTol*(pi/180),
    			dipBandwidth=dipBandwidth)
    		for (i in c(1:length(distanceVector))){
    			pairEquality[i]<-azimuthTest(direction,sampleLocation[pairlist[1,i],],sampleLocation[pairlist[2,i],])
    		}
    	}
    	if (option == 3){
    		if(!(is.null(sampleTime))){
    			timedistVector <- as.vector(vegdist(sampleTime, na.rm=TRUE, method = 'euclidean'))
    			pairEquality <- as.numeric(as.numeric(timedistVector==0) & pairEquality)
    		}
    		else stop("sample dates/times missing")
    	}
    	pairFrame <- cbind(pairEquality,distanceVector, as.vector(simMatrix))
    } else if (option %in% c(2,4)){
    	#Construct pairFrame for locations
		pairEquality <- rep(1,nrow(sampleData)*(nrow(sampleData)-1)*0.5)		
		if (option == 4){
			if (!(is.null(LocationNames))){
    		pairNames <- combn(LocationNames, 2)
    		pairEquality <- apply(pairNames, 2, locationEquality)} 
    	else if(!(is.null(sampleLocation))){
    		pairEquality <- pairEquality & as.vector(vegdist(sampleLocation,na.rm=TRUE,method='euclidean'))
    		pairEquality <- as.numeric(pairEquality==0)
			} else stop("locations missing")	
    	}
    	timeDistanceMatrix <- vegdist(as.numeric(sampleTime), na.rm=TRUE, method='euclidean')
    	timeVector <- as.vector(timeDistanceMatrix)
    	simMatrix <- vegdist(sampleData, method=distmeth, na.rm=TRUE)
	    pairFrame <- cbind(pairEquality, timeVector, as.vector(simMatrix))

    }
	for(i in c(1:(lagNumber+1))){
        pairFrame <- cbind(pairFrame, c(0))
    }
    
    for(j in c(from=1:nrow(pairFrame))) {
    	pairFrame[j,4:dim(pairFrame)[2]] <- distanceLags(lagNumber, lagSize, lagTol, 
    	pairFrame[j,2])
    }
	return(pairFrame)
}

commcorrelogram <- function(sampleData, sampleTime=NULL, sampleLocation=NULL,LocationNames=NULL
	,option = 1,metric='anosim',lagNumber, lagSize, lagTol, numTests=999,anisotropic = FALSE, azimuth
	,azimuthTol,bandwidth,dipAngle,dipTol,dipBandwidth,distmeth='bray'
	,mantmeth='spearman',adj='holm',prog=TRUE,alternative='one.sided') {
#     require(vegan)
	if(!is.null(sampleLocation)){
		if(dim(sampleLocation)[2]<3) stop('Not enough dimensions for xyz coordinates')
	}
  	pairFrame<-lag.sort(sampleData, sampleTime, 		
		sampleLocation,LocationNames,option,lagNumber, 
		lagSize, lagTol, anisotropic, azimuth, 		
		azimuthTol,bandwidth,dipAngle,dipTol,dipBandwidth,distmeth)
    
    #Calculate R-values for lags
    rVals <- RValues(pairFrame,lagNumber,metric,mantmeth)
    
    #Calculate significance
    lagSig <- rep(0,lagNumber+1)
    i<-1
    while(i<=numTests){
        if (!is.null(sampleTime)) sigTime <- sample(sampleTime) else sigTime<-NULL
        if (!is.null(sampleLocation)) 
        	sigLoc <- cbind(sample(sampleLocation[,1])
        	,sample(sampleLocation[,2]),sample(sampleLocation[,3])) else sigLoc<-NULL
        if (!is.null(LocationNames)) sigLN <- sample(LocationNames) else sigLN<-NULL
        sigData <- lag.sort(sampleData, sigTime, 		
			sigLoc,sigLN,option,lagNumber, 
			lagSize, lagTol, anisotropic, azimuth, 		
			azimuthTol,bandwidth,dipAngle,dipTol,dipBandwidth,distmeth)
        testR <- RValues(sigData,lagNumber,metric,mantmeth)
        for(j in c(1:(lagNumber+1))) {
            if (!is.na(rVals[j,1])) {
            	if(alternative=='one.sided'){	
            		if(!(is.na(testR[j,1]))){
            			if(testR[j,1] >= rVals[j,1]) {
                			lagSig[j] <- lagSig[j] + 1
            			}
            		}
				} else if(alternative=='two.sided'){
					if(!(is.na(testR[j,1]))){
            			if(abs(testR[j,1]) >= abs(rVals[j,1])) {
                			lagSig[j] <- lagSig[j] + 1
            			}
            		}
				}
                        } else lagSig[j]<-NA
            
        }
        i<-i+1
    }
    
    lagSig <- lagSig/(numTests+1)
    
    lags<-seq(0,lagNumber*lagSize,lagSize)
    if(prog==F){
    	lagSig<-p.adjust(lagSig,method=adj)
    } else {
    	for (n in c(1:length(lagSig))){
    		temp<-p.adjust(lagSig[1:n],method=adj)
    		lagSig[n]<-temp[n]
    	}
    }
    
	frame<-new('commcorrelogram',community.correlogram=data.frame(rVals[,3],rVals[,2],rVals[,1],lagSig))
    if(adj=='none') names(frame@community.correlogram)<-c("lag.distance","#.Pairs"
    	,paste(metric,"statistic",sep='.'),"significance") else
    names(frame@community.correlogram)<-c("lag.distance","#.Pairs"
    	,paste(metric,"statistic",sep='.'),"adj.significance")
    return(frame)
}


