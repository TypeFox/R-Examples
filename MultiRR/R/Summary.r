median2 <- function(x){apply(x, 2, median)}
lower2 <- function(x){apply(x, 2, quantile, 0.025)}
upper2 <- function(x){apply(x, 2, quantile, 0.975)}
sd2 <- function(x){apply(x, 2, sd)}

Summary <- function(x){
   names2 <- c("SimulatedValue", "Median", "Lower95%", "Upper95%", "Individuals", "Sereies","SeriesPerID", "n.obs")
   n.sim <- x$n.sim[1]
   n.ss <- x$n.ss
   x$SeriesPerInd <- x$Series/x$Individuals

   ##PopulationMeans
   mIntercept <- as.vector(unlist(lapply(x$Intercept, median2)))
   lIntercept <- as.vector(unlist(lapply(x$Intercept, lower2)))
   uIntercept <- as.vector(unlist(lapply(x$Intercept,upper2)))
   Intercept <- data.frame(x$SimInt, mIntercept, lIntercept, uIntercept, x$Individuals, x$Series, x$SeriesPerInd, unlist(x$n.obs))
   colnames(Intercept) <- names2

   mSlope <- as.vector(unlist(lapply(x$Slope, median2)))
   lSlope <- as.vector(unlist(lapply(x$Slope, lower2)))
   uSlope <- as.vector(unlist(lapply(x$Slope,upper2)))
   Slope <- data.frame(x$SimSlope, mSlope, lSlope, uSlope, x$Individuals, x$Series, x$SeriesPerInd, unlist(x$n.obs))
   colnames(Slope) <- names2

   Population.means <- list(Intercept=Intercept, Slope=Slope)

   
   ##Individual VCV
   aID <- lapply(x$IDVCV, median2)
   cID <- lapply(x$IDVCV, lower2)
   bID <- lapply(x$IDVCV,upper2)
  
   medianID <- as.data.frame(matrix(unlist(aID),n.ss,4, byrow=TRUE))
   UpperID  <- as.data.frame(matrix(unlist(bID),n.ss,4, byrow=TRUE))
   LowerID <-  as.data.frame(matrix(unlist(cID),n.ss,4, byrow=TRUE))
   SimVCVInd <- as.data.frame(matrix(unlist(x$SimVCVInd),n.ss,4, byrow=TRUE))
   
    
   IDIntVar <- cbind(SimVCVInd[,1], medianID[,1], LowerID[,1], UpperID[,1], x$Individuals, x$Series, x$SeriesPerInd, unlist(x$n.obs))
   IDCoVar <- cbind(SimVCVInd[,2], medianID[,2], LowerID[,2], UpperID[,2], x$Individuals, x$Series,x$SeriesPerInd, unlist(x$n.obs))
   IDSlopeVar <- cbind(SimVCVInd[,4], medianID[,4], LowerID[,4], UpperID[,4], x$Individuals, x$Series, x$SeriesPerInd, unlist(x$n.obs))
   colnames(IDIntVar) <- names2
   colnames(IDCoVar) <- names2
   colnames(IDSlopeVar) <- names2
   ID <- list(IntVar=IDIntVar, CoVar=IDCoVar, SlopeVar=IDSlopeVar)

   ##Series VCV
   aSeries <- lapply(x$SeriesVCV, median2)
   cSeries <- lapply(x$SeriesVCV, lower2)
   bSeries <- lapply(x$SeriesVCV,upper2)

   medianSeries <- as.data.frame(matrix(unlist(aSeries),n.ss,4, byrow=TRUE))
   UpperSeries  <- as.data.frame(matrix(unlist(bSeries),n.ss,4, byrow=TRUE))
   LowerSeries <-  as.data.frame(matrix(unlist(cSeries),n.ss,4, byrow=TRUE))
   SimVCVSeries <- as.data.frame(matrix(unlist(x$SimVCVSeries),n.ss,4, byrow=TRUE))

   SeriesIntVar <- cbind(SimVCVSeries[,1], medianSeries[,1], LowerSeries[,1], UpperSeries[,1], x$Individuals, x$Series,x$SeriesPerInd, unlist(x$n.obs))
   SeriesCoVar <- cbind(SimVCVSeries[,2], medianSeries[,2], LowerSeries[,2], UpperSeries[,2], x$Individuals, x$Series,x$SeriesPerInd, unlist(x$n.obs))
   SeriesSlopeVar <- cbind(SimVCVSeries[,4], medianSeries[,4], LowerSeries[,4], UpperSeries[,4], x$Individuals,x$Series, x$SeriesPerInd, unlist(x$n.obs))
   colnames(SeriesIntVar) <- names2
   colnames(SeriesCoVar) <- names2
   colnames(SeriesSlopeVar) <- names2
   Series <- list(IntVar=SeriesIntVar, CoVar=SeriesCoVar,SlopeVar=SeriesSlopeVar)


   ##Residuals
   aResiduals <- as.vector(unlist(lapply(x$Residuals, median2)))
   lResiduals <- as.vector(unlist(lapply(x$Residuals, lower2)))
   uResiduals <- as.vector(unlist(lapply(x$Residuals,upper2)))
   Residuals <- as.data.frame(cbind(x$SimResiduals, aResiduals, lResiduals, uResiduals, x$Individuals, x$Series, x$SeriesPerInd, unlist(x$n.obs)))
   colnames(Residuals) <- names2

   
    ##Repeatabilites
    mRInt <- as.vector(unlist(lapply(x$RInt, median)))
    mRSlope <- as.vector(unlist(lapply(x$RSlope, median)))
    mRrn <- as.vector(unlist(lapply(x$Rrn, median)))
    lRInt <- as.vector(unlist(lapply(x$RInt, quantile, 0.025)))
    lRSlope <- as.vector(unlist(lapply(x$RSlope, quantile, 0.025)))
    lRrn <- as.vector(unlist(lapply(x$Rrn, quantile, 0.025)))
    uRInt <- as.vector(unlist(lapply(x$RInt, quantile, 0.975)))
    uRSlope <- as.vector(unlist(lapply(x$RSlope, quantile, 0.975)))
    uRrn <- as.vector(unlist(lapply(x$Rrn, quantile, 0.975)))
   
    SimRepInt <- as.vector(SimVCVInd[,1]/(SimVCVSeries[,1] + SimVCVInd[,1]))
    SimRepSlope <- as.vector(SimVCVInd[,4]/(SimVCVSeries[,4] + SimVCVInd[,4]))
    SimRepRn <- (SimVCVInd[,1] + SimVCVInd[,4] + (2*SimVCVInd[,2]))/((SimVCVSeries[,1] + SimVCVSeries[,4] + (2*SimVCVSeries[,2])) + (SimVCVInd[,1] + SimVCVInd[,4] + (2*SimVCVInd[,2])))
   
    RInt <- as.data.frame(cbind(SimRepInt[1], mRInt, lRInt, uRInt, x$Individuals, x$Series, x$SeriesPerInd, unlist(x$n.obs)))
    colnames(RInt) <- names2

    RSlope <- as.data.frame(cbind(SimRepSlope[1], mRSlope, lRSlope, uRSlope, x$Individuals, x$Series, x$SeriesPerInd, unlist(x$n.obs)))
    colnames(RSlope) <- names2
    
    Rrn <- as.data.frame(cbind(SimRepRn[1], mRrn, lRrn, uRrn, x$Individuals, x$Series, x$SeriesPerInd, unlist(x$n.obs)))
    colnames(Rrn) <- names2
    Repeatabilities <- list(Intercept=RInt, Slope=RSlope)
                     
   summary <- list(Population.means=Population.means, ID=ID, Series=Series, Repeatabilities=Repeatabilities, Residuals=Residuals)
   print(summary)
}
