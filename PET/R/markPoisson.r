#########################################################################
#
# Copyright Weierstrass Institute for Applied Analysis and 
#           Stochastics (WIAS) & Humboldt Universitaet zu Berlin, 
#           Germany 2006
# *********************************************************
#
# Name:          markPoisson.m
#                ---------------
# Author:        Joern Schulz
# Stand:         13.08.2006
#
# Uses:          
#
#########################################################################

markPoisson <- function(DataInt, nSample=200000, 
				ThetaSamples=181,
				RhoSamples=2*round(sqrt(sum((dim(DataInt))^2))/2)+1,
				RhoMin=-0.5*sqrt(2),
                        vect.length = 100000, image=TRUE, DebugLevel="Normal"){

# old default: 
#    RhoSamples=2*ceiling(sqrt(sum((dim(DataInt)-floor((dim(DataInt)-1)/2)-1)^2)))+3,
#######################################################################
#
# DataInt   (matrix)  An two dimensional image in which the matrix elements 
#           comply to different intensities. That mean, how intense is the 
#           decay of positrons in a certain tissue.
# nSample   (integer)  nSample determine the number of accepted events which 
#           to be generated with AR method. Defaults to nSample = 200000.
# ThetaSamples (integer)  Specifies the number of samples in the angular 
#           parameter theta in the sinogram. The sinogram is sampled linearly 
#           from 0 to (approximately) \eqn{\pi}{pi} radians. Defaults to
#           ThetaSamples=181.
# RhoSamples (integer) Specifies the number of samples in the distance 
#           parameter rho in the sinogram. Defaults to RhoSamples =
#           2*ceiling(sqrt(sum((dim(oData) - floor((dim(DataInt)-1)/2)-1)^2)))+3.
# vect.length (integer) Determine a bound of number of generated accepted 
#           events in each iteration. That mean, if nSample>vect.length then 
#           in each iteration will be generate vect.length/(0.1217) events 
#           and roughly vect.length accepted events. If you scale down 
#           vect.length then the computation is more time-consuming and if 
#           you increase value of vect.length then more memory is necessary.
# image     (logical)  If image=FALSE then only an generated sinogram returns, 
#           i.e. a noisy Radon Transform of the phantom returns. Otherwise, if 
#           image=TRUE additionaly an image returns, where each pixel of the 
#           image corresponds to number of accepted events generated with the 
#           AR method. Defaults to \code{image=TRUE}. }
# DebugLevel (character) This parameter controls the level of output. Following
#           possibilities are available: The default "Normal" for standard level 
#           of output to screen or alternative "Detail" if it desirable to logged 
#           allmost all output to screen or "HardCore" for no information at all.
#
#######################################################################
#
#
  args <- match.call()
  
  pcType <- .Platform$OS.type # discern of '"unix"' or '"windows"'
  if (pcType=="unix") pcType <- 1
  else if (pcType=="windows") pcType <- 2
  else pcType <- 0
  
  DL1 <- logDebug(DebugLevel)
  DebugLevel <- DL1[[3]]
  DL2 <- DL1[[2]]
  DL1 <- DL1[[1]]
  
  RhoMin <- abs(RhoMin)                     
  M <- nrow(DataInt)
  N <- ncol(DataInt)
  T <- ThetaSamples
  R <- RhoSamples 

  ptmAll <- proc.time()

  maxDataInt <- max(DataInt)		# constant C for AR-method is specify
  PData  <- matrix(0, nrow=M , ncol=N )   # initialization of PData
  rPData <- matrix(0, nrow=T , ncol=R )   # initialization of rPData

  if(DL1) cat("AR step: \n")

  if(nSample > vect.length ){
    nRemainder <- nSample
    nSample <- vect.length

    if(DL2) cat("To generating events: ", nRemainder,"\n")

    repeat{
      n <- ceiling(nSample/(0.1217))	# n is the expectation value, which is necessary for the computation of nSample-variables. 0.1217 was calculated from different tests.

      if(DL2) cat("AR step: \n")
      ptm <- proc.time()

     # ACCEPTANCE-REAJECTION METHOD
     # Generating of random events with specific properties
      xCord <- runif(n)				# Generating of X
      yCord <- runif(n)				# Generating of Y
      Samp <- matrix(c(xCord,yCord),ncol=2) # Collecting in a (n,2)-matrix

     # -------> ACCEPTANCE-REAJECTION-STEP <-------
      U <- runif(n)				# Generating of U
  
     # Assignment of the events to voxels:
      SampVox <- matrix(c(ceiling(M*xCord),ceiling(N*yCord)),ncol=2) 
      #ni <- 1:n
      #V1indexInt <- na.omit( ifelse( DataInt[index[ni,]]/maxDataInt > U, ni, NA ))
      SampAcceptR <- DataInt[SampVox]/maxDataInt > U	# acceptance?
      SampAccept <- Samp[SampAcceptR,]			# Selection of permissible random variables
    
     # if too much random variables was generated then choose the first 'nRemainder'
      if(nrow(SampAccept) > nRemainder) SampAccept <- SampAccept[1:nRemainder,]
      if(DL2){
        cat("Needed time for creation of accepted events: ", proc.time()-ptm, "\n")
        cat("Number of generated events: ", nrow(SampAccept), "\n")
      }

      ######################################################

     # Remove of not needed variables
      rm(xCord, yCord, Samp, SampVox, SampAcceptR) 

      if(image){
       # Summarization of the generated events in a (M,N)-matrix. 
       # The summarization happens through a classification of [0,1] in M-intervals for
       # X1 and N-intervals for X2
       #ptm <- proc.time()				# question time
        xSampVox <- cut(SampAccept[,1], seq(0,1,length.out=M+1))
        ySampVox <- cut(SampAccept[,2], seq(0,1,length.out=N+1))
        PData <- PData + table(xSampVox,ySampVox) # Summarization of the events
        rm(xSampVox,ySampVox)
        #cat("Benoetigte Zeit zur  Ereignisszuweisung: ", proc.time()-ptm, "\n")
      }

      if(DL2) cat("Computation of parameters for Radon Transformation: \n")
      ptm <- proc.time()
      lSampA <- nrow(SampAccept)
     # Generating of  angles to the accepted events,
      theta <- runif(lSampA, min=0, max=pi)

     # computation the distances
     # Scaling the image for the computation of the distances, so that X,Y is in
     # [-0.5,0.5]. The centre of image is than [0,0].
      SampAccept <- SampAccept-0.5
      radonSamp <- SampAccept[,1]*cos(theta)+SampAccept[,2]*sin(theta)

     # Subdivision in intervals
      radonSamp <- cut(radonSamp, seq(-RhoMin,RhoMin,length.out=R+1))
      theta <- cut(theta, pi/T*(0:T))
      #radonSamp <- cbind(theta,radonSamp)

      if(DL2) cat("Needed time for creation of theta and computation of rho: ", proc.time()-ptm, "\n")

     # Cumulation of equal intervals 
      #rPData <- table(radonSamp[,1], radonSamp[,2])
      rPData <- rPData + table(theta, radonSamp)

      nRemainder <- nRemainder - lSampA
      
      if(DL1 & !DL2){
        if (pcType==1){
          cat("Progress:",trunc(sum(rPData)*100/(sum(rPData)+nRemainder)),"%\r");
        } else if (pcType==2) {
          cat("Progress:",trunc(sum(rPData)*100/(sum(rPData)+nRemainder)),"%\r");
          flush.console();
        }
      }
      if(DL2)
        cat("Remaining number of to generating events:", nRemainder,"\n")
    
      if( nRemainder <= vect.length) break
      nSample <- vect.length
    } 

  nSample <- nRemainder
  }


######################################################################################
######################################################################################

  if(nSample!=0){
    SampAccept <- NULL
    n <- ceiling(nSample/(0.1217))	# n is the expectation value, which is necessary for the computation of nSample-variables. 0.1217 was calculated from different tests.

    if(DL2) cat("AR step: \n")
    ptm <- proc.time()

    repeat{
     # ACCEPTANCE-REAJECTION METHOD
     # Generating of random events with specific properties
      xCord <- runif(n)				# Generating von X
      yCord <- runif(n)				# Generating von Y
      Samp <- matrix(c(xCord,yCord),ncol=2) # Collecting in a (n,2)-Matrix

     # -------> ACCEPTANCE-REAJECTION-STEP <-------
      U <- runif(n)				# Generating of U
  
     # Assignment of the events to voxels:
      SampVox <- matrix(c(ceiling(M*xCord),ceiling(N*yCord)),ncol=2) 
      #ni <- 1:n
      #V1indexInt <- na.omit( ifelse( DataInt[index[ni,]]/maxDataInt > U, ni, NA ))
      SampAcceptR <- DataInt[SampVox]/maxDataInt > U	# acceptance?
      SampAcceptR <- Samp[SampAcceptR,]			# Selection of permissible random variables
      SampAccept <- rbind(SampAccept,SampAcceptR)
  
      if(nrow(SampAccept) >= nSample) break
      n <- ceiling((nSample-nrow(SampAccept))/(0.1217))
    }
   # if too much random variables was generated then choose the first 'nSample'
    if(nrow(SampAccept) > nSample) SampAccept <- SampAccept[1:nSample,]
    
    if(DL2){
      cat("Needed time for creation of accepted events:  ", proc.time()-ptm, "\n")
      cat("Number of generated events: ", nrow(SampAccept), "\n")
    }

  ######################################################

    # Remove of not needed variables
    rm(xCord, yCord, Samp, SampVox, SampAcceptR) 

    if(image){
     # Summarization of the generated events in a (M,N)-matrix. 
     # The summarization happens through a classification of [0,1] in M-intervals for
     # X1 and N-intervals for X2
     #ptm <- proc.time()					# question time
      xSampVox <- cut(SampAccept[,1], seq(0,1,length.out=M+1))
      ySampVox <- cut(SampAccept[,2], seq(0,1,length.out=N+1))
      PData <- PData + table(xSampVox,ySampVox)		# Summarization of the events
      rm(xSampVox,ySampVox)
      #cat("Benoetigte Zeit zur  Ereignisszuweisung: ", proc.time()-ptm, "\n")
    }

    if(DL2) cat("Computation of parameters for Radon Transformation: \n")
    ptm <- proc.time()
    lSampA <- nrow(SampAccept)
   # Generating of  angles to the accepted events.
    theta <- runif(lSampA, min=0, max=pi)

   # computation the distances
   # Scaling the image for the computation of the distances, so that X,Y is in
   # [-0.5,0.5]. The centre of image is than [0,0].
    SampAccept <- SampAccept-0.5
    radonSamp <- SampAccept[,1]*cos(theta)+SampAccept[,2]*sin(theta)

   # Subdivision in intervals
    radonSamp <- cut(radonSamp, seq(-RhoMin,RhoMin,length.out=R+1))
    theta <- cut(theta, pi/T*(0:T))
    #radonSamp <- cbind(theta,radonSamp)

    if(DL2) cat("Needed time for creation of theta and computation of rho: ", proc.time()-ptm, "\n")

   # Cumulation of equal intervals
    #rPData <- table(radonSamp[,1], radonSamp[,2])
    rPData <- rPData + table(theta, radonSamp)
    
    if(DL1 & !DL2){
      if (pcType==1){
        cat("Progress:",trunc(100),"%\r");
      } else if (pcType==2) {
        cat("Progress:",trunc(100),"%\r");
        flush.console();
      }
    }
    }

######################################################################################
######################################################################################

  if(DL1)  
    cat("Needed time in all: ", proc.time()-ptmAll, "\n")

  #DeltaXY <- c(pi/ThetaSamples, 1/sqrt(2))
  #XYmin <- c(0, -0.5*(DeltaXY[2]*(RhoSamples-1)))
  XYmin <- c(0, (-0.5*((2*round(sqrt(sum((dim(DataInt))^2))/2)+1)-1)))
  DeltaXY <- c(pi/ThetaSamples, (2*abs(XYmin[2])+1)/RhoSamples)
  if(image) {
      z <- list(rData=rPData, Data=PData,
                Header=list(SignalDim=c(ThetaSamples,RhoSamples), 
                           XYmin=XYmin, 
                           DeltaXY=DeltaXY), 
                call=args)
  } else { 
      z <- list(rData=rPData,
                Header=list(SignalDim=c(ThetaSamples,RhoSamples), 
                            XYmin=XYmin, 
                            DeltaXY=DeltaXY), 
                call=args)
  }
  class(z)<-"pet"
  z
  }