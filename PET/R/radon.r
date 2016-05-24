#########################################################################
#
# Copyright Weierstrass Institute for Applied Analysis and 
#           Stochastics (WIAS) & Humboldt Universitaet zu Berlin, 
#           Germany 2006
# *********************************************************
#
# Name:          PETphantom.m
#                ---------------
# Author:        Joern Schulz
# Stand:         07.08.2006
#
#########################################################################

radon <- function(oData, 
                  mode="NN", 
			XYSamples=nrow(oData),
			XYmin=-0.5*(nrow(oData)-1),
			DeltaXY=1,
			ThetaSamples=181,
			RhoSamples=2*round(sqrt(sum((dim(oData))^2))/2)+1,
			ThetaMin=0,
                  RhoMin=-0.5*((2*round(sqrt(sum((dim(oData))^2))/2)+1)-1),
			DeltaTheta=pi/ThetaSamples,
			DeltaRho=(2*abs(RhoMin)+1)/RhoSamples)
{
# old default: 
#    RhoSamples=2*ceiling(sqrt(sum((dim(oData)-floor((dim(oData)-1)/2)-1)^2)))+3,
#    RhoMin=-0.5*(RhoSamples-1),
#    DeltaRho=1
#    ThetaSamples=ceiling(pi*(nrow(oData)-1))
#    RhoMin=-0.5*(RhoSamples-1)*sqrt(2),
#    DeltaRho=1*sqrt(2)
#

# ========================================================================
#
# XYmin           Specifies the minimum sample position in the image on the 
#                first axis. If nor given the image is centered around the 
#                middle.
# Ymin           Specifies the minimum sample position in the image on the 
#                second axis. If nor given the image is centered around the 
#                middle.
# DeltaXY         Specifies the sampling distance of both axes in the image.
# ThetaSamples   Specifies the number of samples in the angular parameter
#                theta in the sinogram. The sinogram is sampled linearly
#                from 0 to (approximately) pi radians. The default sampling 
#                distance in theta is pi/ThetaSamples
# RhoSamples     Specifies the number of samples in the distance parameter
#                rho in the sinogram. Should be an odd number.
# DeltaRho       Specifies the sampling distance in rho. The program will center
#                the sampling points around 0.
# ThetaMin       Specifies the minimum sample position in the sinogram on 
#                the first axis.
# RhoMin         Specifies the minimum sample position in the sinogram on 
#                the second axis.
#
# Examples:
#
#  A <- phantom()
#  R <- radon(A$data)
#  image(R$rdata,col=gray((0:255)/255))
#  R <- radon(A$data, DeltaRho=1/sqrt(4))
#  image(R$rdata,col=gray((0:255)/255))
#

  args <- match.call()
 
 ######################################### 
 # checking parameter
  #if (DeltaRho>(DeltaXY/sqrt(2)))
  #    stop("DeltaRho should be less than DeltaXY/sqrt(2)")
  if (!(is.matrix(oData)))
      stop("'oData' has to be of type 'matrix'.")
  if (nrow(oData)!=ncol(oData))
      stop("At the moment only image of quadratic form are admissible.")

 ######################################### 
 # setting all parameters to setpar

 # Fundamental setting.
  setpar      <- matrix(0, nrow=9, ncol=1)
  setpar[1,1] <- XYSamples
  setpar[2,1] <- XYmin
  setpar[3,1] <- DeltaXY
  setpar[4,1] <- ThetaSamples 
  setpar[5,1] <- DeltaTheta
  setpar[6,1] <- RhoSamples
  setpar[7,1] <- RhoMin
  setpar[8,1] <- DeltaRho
  setpar[9,1] <- ThetaMin


 ###########################################################
 # calling the C-Routines
  if (mode=="NN") {
    
      rdata <- .C("radonNN", 
                  rdata=matrix(0, nrow=ThetaSamples, ncol=RhoSamples), 
                  as.double(oData), 
                  as.double(setpar),
                  PACKAGE="PET")$rdata

  } else if(mode=="LI"){

      rdata <- .C("radonLI", 
                  rdata=matrix(0, nrow=ThetaSamples, ncol=RhoSamples), 
                  as.double(oData), 
                  as.double(setpar),
                  PACKAGE="PET")$rdata

  } else if(mode=="SINC"){

      rdata <- .C("radonSINC", 
                  rdata=matrix(0, nrow=ThetaSamples, ncol=RhoSamples), 
                  as.double(oData), 
                  as.double(setpar),
                  PACKAGE="PET")$rdata

  } else
      stop("The mode=",mode," is not supported.")

  

  z <- list(rData=rdata, 
            Header=list(SignalDim=c(ThetaSamples, RhoSamples), 
                        XYmin=c(ThetaMin, RhoMin), 
                        DeltaXY=c(DeltaTheta, DeltaRho)), 
            call=args )
  class(z) <- "pet"
  return(z)
}
