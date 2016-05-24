Bias <-
function(x){
   names2 <- c("Vint", "COVintslope", "COVslopeint", "Vslope")
   n.sim <- x$n.sim[1]
   n.ss <- x$n.ss

   aID <- lapply(x$IDVCV, median2)
   
   medianID <- as.data.frame(matrix(unlist(aID),n.ss,4, byrow=TRUE))
   SimVCVInd <- as.data.frame(matrix(unlist(x$SimVCVInd),n.ss,4, byrow=TRUE))
   BiasID <- as.data.frame(medianID-SimVCVInd)
   RelBiasID <- as.data.frame(abs(medianID-SimVCVInd)/abs(SimVCVInd))*100
   colnames(medianID) <- names2
   medianID$Individuals <- x$Individuals
   medianID$SeriesPerID <- x$SeriesPerIndividual
   colnames(RelBiasID) <- names2
   RelBiasID$Individuals <- x$Individuals
   RelBiasID$Series <- x$Series
   RelBiasID$SeriesPerID <- x$Series/x$Individuals
   RelBiasID$n.obs <-unlist(x$n.obs)
   colnames(BiasID) <- names2
   BiasID$Individuals <- x$Individuals
   BiasID$Series <- x$Series
   BiasID$SeriesPerID <- x$Series/x$Individuals
   BiasID$n.obs <-unlist(x$n.obs)

####
   aSeries <- lapply(x$SeriesVCV, median2)
   medianSeries <- as.data.frame(matrix(unlist(aSeries),n.ss,4, byrow=TRUE))
   SimVCVSeries <- as.data.frame(matrix(unlist(x$SimVCVSeries),n.ss,4, byrow=TRUE))
   BiasSeries <- as.data.frame(medianSeries-SimVCVSeries)
   RelBiasSeries <- as.data.frame(abs(medianSeries-SimVCVSeries)/abs(SimVCVSeries))*100
   colnames(medianSeries) <- names2
   medianSeries$Individuals <- x$Individuals
   medianSeries$SeriesPerID <- x$SeriesPerIndividual
   colnames(RelBiasSeries) <- names2
   RelBiasSeries$Individuals <- x$Individuals
   RelBiasSeries$Series <- x$Series
   RelBiasSeries$SeriesPerID <- x$Series/x$Individuals
   RelBiasSeries$n.obs <-unlist(x$n.obs)
   colnames(BiasSeries) <- names2
   BiasSeries$Individuals <- x$Individuals
   BiasSeries$Series <- x$Series
   BiasSeries$SeriesPerID <- x$Series/x$Individuals
   BiasSeries$n.obs <-unlist(x$n.obs)
####
RInt <- lapply(x$RInt, median)
RSlope <- lapply(x$RSlope, median)

RepInt <- as.data.frame(matrix(unlist(RInt),n.ss,1, byrow=TRUE))
RepSlope <- as.data.frame(matrix(unlist(RSlope),n.ss,1, byrow=TRUE))

   
SimRepInt <- SimVCVInd[,1]/(SimVCVSeries[,1] + SimVCVInd[,1])
SimRepSlope <- SimVCVInd[,4]/(SimVCVSeries[,4] + SimVCVInd[,4])
SimRepRn <- (SimVCVInd[,1] + SimVCVInd[,4] + (2*SimVCVInd[,2]))/((SimVCVSeries[,1] + SimVCVSeries[,4] + (2*SimVCVSeries[,2])) + (SimVCVInd[,1] + SimVCVInd[,4] + (2*SimVCVInd[,2])))
   
BiasRepInt1 <- (RepInt-SimRepInt)
RelBiasRepInt1 <- (abs(RepInt-SimRepInt)/SimRepInt)*100

BiasRepSlope1 <- (RepSlope-SimRepSlope)
RelBiasRepSlope1 <- (abs(RepSlope-SimRepSlope)/SimRepSlope)*100

   
BiasRepInt <- data.frame(Repeatability=BiasRepInt1, Individuals=x$Individuals,  Series=x$Series, SeriesPerIndividual=x$Series/x$Individuals, n.obs=unlist(x$n.obs))
colnames(BiasRepInt)[1] <-"Repeatabilty" 

BiasRepSlope <- data.frame(Repeatability=BiasRepSlope1, Individuals=x$Individuals,  Series=x$Series, SeriesPerIndividual=x$Series/x$Individuals, n.obs=unlist(x$n.obs)) 
colnames(BiasRepSlope)[1] <-"Repeatabilty" 

     
RelBiasRepInt <- data.frame(Repeatability=RelBiasRepInt1, Individuals=x$Individuals,  Series=x$Series, SeriesPerIndividual=(x$Series/x$Individuals), n.obs=unlist(x$n.obs))
colnames(RelBiasRepInt)[1] <-"Repeatabilty"    

RelBiasRepSlope <- data.frame(Repeatability=RelBiasRepSlope1, Individuals=x$Individuals,  Series=x$Series, SeriesPerIndividual=x$Series/x$Individuals, n.obs=unlist(x$n.obs)) 
colnames(RelBiasRepSlope)[1] <-"Repeatabilty"    



RepInt <- list(Bias=round(BiasRepInt,2), RelBias=round(RelBiasRepInt,2))
RepSlope <- list(Bias=round(BiasRepSlope,2), RelBias=round(RelBiasRepSlope,2))

Rep <- list(Int=RepInt, Slope=RepSlope)
   
###

   ID <- list( Bias=round(BiasID,2), RelBias=round(RelBiasID,2))
   Series <- list(Bias=round(BiasSeries,2),RelBias=round(RelBiasSeries,2) )
   Bias <- list(SimID=x$SimVCVInd[[1]], SimSeries=x$SimVCVSeries[[1]], SimRepInt=SimRepInt[1],SimRepSlope=SimRepSlope[1], ID=ID, Series=Series, Rep=Rep)
   Bias
   print(Bias)
}
