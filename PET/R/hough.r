hough <- function(oData,
                  mode=1,
                  XYSamples=nrow(oData),
                  DeltaXY=1.0,
                  XYmin=-0.5*DeltaXY*(XYSamples-1),
                  ThetaSamples=181,
                  RhoSamples=2*round(sqrt(sum((dim(oData))^2))/2)+1,
                  ThetaMin=0,
                  RhoMin=-0.5*((2*round(sqrt(sum((dim(oData))^2))/2)+1)-1),
                  DeltaTheta=pi/ThetaSamples,
                  DeltaRho=(2*abs(RhoMin)+1)/RhoSamples,
                  DebugLevel="Normal")
{
# old default: 
#    RhoSamples=2*ceiling(sqrt(sum((dim(oData)-floor((dim(oData)-1)/2)-1)^2)))+3,
#    DeltaRho=1.0,
#    RhoMin=-0.5*(DeltaRho*(RhoSamples-1)),
 # =====================================================================
 # check and setting input parameter

  DL1 <- logDebug(DebugLevel)
  DebugLevel <- DL1[[3]]
  DL2 <- DL1[[2]]
  DL1 <- DL1[[1]]
    
  #if (DeltaRho>(DeltaXY/sqrt(2)))
  #    stop("DeltaRho should be less than DeltaXY/sqrt(2)")
  if (!(is.matrix(oData)))
      stop("'oData' has to be of type 'matrix'.")
  if (nrow(oData)!=ncol(oData))
      stop("At the moment only image of quadratic form are admissible.")

 # =====================================================================
 # call the c-routine

 # Fundamental setting.
  setpar      <- matrix(0, nrow=9, ncol=1)
  setpar[1,1] <-   XYSamples     
  setpar[2,1] <-   XYmin
  setpar[3,1] <-   DeltaXY
  setpar[4,1] <-   ThetaSamples 
  setpar[5,1] <-   DeltaTheta
  setpar[6,1] <-   RhoSamples
  setpar[7,1] <-   RhoMin
  setpar[8,1] <-   DeltaRho
  setpar[9,1] <-   ThetaMin


  if (mode == 1){
     if (DL1) cat("Using Hough1 --> ")
     rdata <- .C("Hough1", 
                 y = matrix(0, nrow=setpar[4,1], ncol=setpar[6,1]), 
                 as.double(oData), 
                 as.double(setpar), 
                 PACKAGE="PET")$y
  } else if (mode == 2){
     if (DL1) cat("Using Hough2 --> ")
     rdata <- .C("Hough2", 
                 y = matrix(0, nrow=setpar[4,1], ncol=setpar[6,1]), 
                 as.double(oData), 
                 as.double(setpar),
                 PACKAGE="PET")$y
  } else if (mode == 3){
     if (DL1) cat("Using Hough3 --> ")
     rdata <- .C("Hough3", 
                 y = matrix(0, nrow=setpar[4,1], ncol=setpar[6,1]), 
                 as.double(oData), 
                 as.double(setpar),
                 PACKAGE="PET")$y
  } else if (mode == 4){
     if (DL1) cat("Using Hough4 --> ")
     rdata <- .C("Hough4", 
                 y = matrix(0, nrow=setpar[4,1], ncol=setpar[6,1]), 
                 as.double(oData), 
                 as.double(setpar),
                 PACKAGE="PET")$y
  } else stop("Only option 'mode' = 1, 2, 3 or 4 possible.")
  
  if (DL1) cat("complete \n")
  
  z <- list(hData=rdata, 
            Header=list(SignalDim=c(ThetaSamples, RhoSamples), 
                        XYmin=c(ThetaMin, RhoMin), 
                        DeltaXY=c(DeltaTheta, DeltaRho)),
            call=args)
  class(z) <- "pet"
  return(z)
}
