ParetoDensityEstimation= function(Data,paretoRadius=NULL,kernels=NULL,MinAnzKernels=100){
#  V = ParetoDensityEstimation(Data,ParetoRadius,Kernels)
#  V= ParetoDensityEstimation(Data,ParetoRadius)
#  Estimates the Pareto Density for a one dimensional distibution
#  this is the best density estimation to judge Gaussian Mixtures  of the Data see [Ultsch 2003}
# 
#  INPUT
#  Data                    die eindimensional verteilten Daten
#  OPTIONAL
#  paretoRadius            der Pareto Radius, wenn nicht angegeben, wird er berechnet
#  kernels                 Data values at which ParetoDensity is measured , use plot(Kernels,ParetoDensity) for display
#                          wird bestimmt, wenn nicht angegeben oder Kernels ==0
#  MinAnzKernels           Minimale Anzahl Kernls, wenn nicht angegeben oder MinAnzKernelss ==0 =>  MinAnzKernels==100	
# 
#  OUTPUT
#  Kernels                 Data values at which ParetoDensity is measured , use plot(Kernels,ParetoDensity) for display
#  paretoDensity           die mit dem ParatoRadius ermittelte Dichte
#  paretoRadius            der Pareto Radius
# 
#  Author: MT 2015/06 uebernommen von ALU 2003

#require(caTools)
# function noNaN_hlp
###############################################
noNaN_hlp<-function(data){
  # Data=noNaN_hlp(data)
  # Entfernt NANs und gibt auch die Indizies wieder
  #
  # INPUT
  # x[,d]     Datensatz der Spaltendimension d
  #
  # OUTPUT
  # elements      
  # noNaNInd
  # nrElements

  # Autor: ALU
  # 1. Editor: MT
  variables <- ncol(data)	
if(!length(variables)){ #MT: Falls statt Matrix Vektor uebergeben wird 
	noNaNInd<- which(is.finite(data))
	elements<-data[noNaNInd]
	nrElements<-length(noNaNInd)
}else{
  noNaNIndAll=c()
	for(i in 1:variables){
		noNaNInd<- which(is.finite(data[,i]))
		data<-data[noNaNInd,]	
		noNaNIndAll=cbind(noNaNIndAll,noNaNInd)
		noNaNInd=c()
	}  
	elements=data
	nrElements<-nrow(noNaNIndAll)
  noNaNInd=noNaNIndAll
}
return(list(elements=elements,noNaNInd=noNaNInd,nrElements=nrElements))
}
###############################################



	if(paretoRadius==0 || length(paretoRadius)==0)
		paretoRadius=ParetoRadius(Data)
	
	Data <- noNaN_hlp(Data)$elements
	Data <- as.matrix(Data)


	r <- nrow(Data)
	c <- ncol(Data)



	if(c>r)
		Data <- t(Data)


	if(c>1)
		return('ERROR ParetoDensityEstimation: Data set not univariate !')

	if(length(kernels)==0 || kernels==0){ #MT: Korrektur,?: statt kernels==0 und im Input Kernels=0
		nBins <-OptimalNoBins(Data[,1])
		#MT: MinAnzKernels fehlte
		nBins = max(MinAnzKernels ,nBins);  # mindestzahl von Kernels sicherstellen
		if (nBins>100){ 
		  if(nBins>1E4){#MT: Fehlerabdfang bei zu vielen Bins
		    rHist <- hist(Data,seq(min(Data), max(Data), length.out=1E4),plot=FALSE)
		    warning('Too many bins estimated, try to transform or sample the data')
		  }else{
			  rHist <- hist(Data,seq(min(Data), max(Data), length.out=(nBins*3)+1),plot=FALSE)
		  }
		}else{ #MT: Dadurch ist bei <100 die Formel nBins*3+1 unn?tig, da nBins=100 gesetzt
			rHist <- hist(Data,seq(min(Data), max(Data), length.out=nBins),plot=FALSE)
		}
		nInBins <- rHist$counts
		kernels <- rHist$mids
	}
	
	nKernels <- length(kernels)

	minData <- min(Data,na.rm=TRUE)
	maxData <- max(Data,na.rm=TRUE)

#  diese Daten liegen am unteren Rand
	lowBInd <-  (Data < (minData+paretoRadius) )
	lowR <- as.matrix(2*minData-Data[lowBInd],ncol=1) #RG: Data[lowBInd,] falsche Dimension, da Data eindimensional
	
# diese Daten liegen am obere Rand
	upBInd <-  (Data > (maxData-paretoRadius) )
	upR <- as.matrix(2*maxData-Data[upBInd],ncol=1)


	DataPlus <- as.matrix(c(Data,lowR,upR),1)
	
	paretoDensity <- rep(0,nKernels)
	
	for(i in 1:nKernels){
		lb <- kernels[i] - paretoRadius
		ub <- kernels[i] + paretoRadius
		isInParetoSphere <- (DataPlus >= lb) & (DataPlus <= ub)
		paretoDensity[i] <- sum(isInParetoSphere)
	}
	
# trapz funktion (caTools)
	area <- trapz(kernels,paretoDensity)
	#idx = 2:length(kernels)
    #area <- (as.double((kernels[idx] - kernels[idx - 1]) %*% (paretoDensity[idx] + paretoDensity[idx - 1]))/2)
	
	if(area < 0.0000000001 || is.na(area)) #RG: Fall: kernel==0 => area==NAN muss abgefangen werden
		paretoDensity <- rep(0,nKernels)
	else 
		paretoDensity <- paretoDensity /area

 return (list(kernels=kernels,paretoDensity=paretoDensity,paretoRadius=paretoRadius)) 

 }

