Imprecision <-
function(x){
   names2 <- c("Vint", "COVintslope", "COVslopeint", "Vslope")
   n.sim <- x$n.sim[1]
   n.ss <- x$n.ss
       
   aID <- lapply(x$IDVCV, mean2)
   bID <- lapply(x$IDVCV, sd2)
   
   meanID <- as.data.frame(matrix(unlist(aID),n.ss,4, byrow=TRUE))
   sdID <- as.data.frame(matrix(unlist(bID),n.ss,4, byrow=TRUE))
   ImprecisionID <- as.data.frame(abs(sdID/meanID))*100

   colnames(ImprecisionID) <- names2
   ImprecisionID$Individuals <- x$Individuals
   ImprecisionID$Series <- x$Series
   ImprecisionID$SeriesPerID <- x$Series/x$Individuals
   ImprecisionID$n.obs <- unlist(x$n.obs)

    SimVCVInd <- as.data.frame(matrix(unlist(x$SimVCVInd),n.ss,4, byrow=TRUE))
    SimVCVSeries <- as.data.frame(matrix(unlist(x$SimVCVSeries),n.ss,4, byrow=TRUE))
####
   aSeries <- lapply(x$SeriesVCV, mean2)
   bSeries <- lapply(x$SeriesVCV, sd2)
   
   meanSeries <- as.data.frame(matrix(unlist(aSeries),n.ss,4, byrow=TRUE))
   sdSeries <- as.data.frame(matrix(unlist(bSeries),n.ss,4, byrow=TRUE))
   ImprecisionSeries <- as.data.frame(abs(sdSeries/meanSeries))*100

   colnames(ImprecisionSeries) <- names2
   ImprecisionSeries$Individuals <- x$Individuals
   ImprecisionSeries$Series <- x$Series
   ImprecisionSeries$SeriesPerID <- x$Series /x$Individual

   ImprecisionSeries$n.obs <- unlist(x$n.obs)

###
   SimRepInt <- SimVCVInd[,1]/(SimVCVSeries[,1] + SimVCVInd[,1])
   SimRepSlope <- SimVCVInd[,4]/(SimVCVSeries[,4] + SimVCVInd[,4])
   SimRepRn <- (SimVCVInd[,1] + SimVCVInd[,4] + (2*SimVCVInd[,2]))/((SimVCVSeries[,1] + SimVCVSeries[,4] + (2*SimVCVSeries[,2])) + (SimVCVInd[,1] + SimVCVInd[,4] + (2*SimVCVInd[,2])))
   
   mRInt <- lapply(x$RInt, mean)
   mRSlope <- lapply(x$RSlope, mean)
   mRrn <- lapply(x$Rrn, mean)
   
   sdRInt <- lapply(x$RInt, sd)
   sdRSlope <- lapply(x$RSlope, sd)
   sdRrn <- lapply(x$Rrn, sd)
   
   
   meanRInt <- as.data.frame(matrix(unlist(mRInt),n.ss,1, byrow=TRUE))
   sdRInt <- as.data.frame(matrix(unlist(sdRInt),n.ss,1, byrow=TRUE))
   ImprecisionRInt <- as.data.frame(sdRInt/meanRInt)*100

   colnames(ImprecisionRInt) <- "Repeatability"
   ImprecisionRInt$Individuals <- x$Individuals
   ImprecisionRInt$Series <- x$Series
   ImprecisionRInt$SeriesPerID <- x$Series/x$Individual
   ImprecisionRInt$n.obs <- unlist(x$n.obs)
   
   meanRSlope <- as.data.frame(matrix(unlist(mRSlope),n.ss,1, byrow=TRUE))
   sdRSlope <- as.data.frame(matrix(unlist(sdRSlope),n.ss,1, byrow=TRUE))
   ImprecisionRSlope <- as.data.frame(sdRSlope/meanRSlope)*100

   colnames(ImprecisionRSlope) <- "Repeatability"
   ImprecisionRSlope$Individuals <- x$Individuals
   ImprecisionRSlope$Series <- x$Series
   ImprecisionRSlope$SeriesPerID <- x$Series/x$Individual
   ImprecisionRSlope$n.obs <- unlist(x$n.obs)


   
Imprecision <- list(SimID=x$SimVCVInd[[1]], SimSeries=x$SimVCVSeries[[1]], SimRepInt=SimRepInt, SimRepInt=SimRepRn, SimRepSlope=SimRepSlope[1], SimRepInt=SimRepRn, ImprecisionID=ImprecisionID, ImprecisionSeries=ImprecisionSeries, ImprecisionRInt=ImprecisionRInt, ImprecisionRSlope=ImprecisionRSlope)
print(Imprecision)
}
