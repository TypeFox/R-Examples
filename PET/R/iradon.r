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
# Stand:         30.03.2006
#
# Uses:          
#
#########################################################################
# 
#
iradon <- function(rData, XSamples, YSamples, 
                   mode="FB", 
                   Xmin=-sqrt(0.5)*(ncol(rData)/XSamples)*0.5*(XSamples-1), 
                   Ymin=-sqrt(0.5)*(ncol(rData)/YSamples)*0.5*(YSamples-1), 
                   DeltaX=sqrt(0.5)*(ncol(rData)/XSamples), 
                   DeltaY=sqrt(0.5)*(ncol(rData)/YSamples),
                   InterPol=1,
                   FilterTyp="Hamming1",
                   oData=NULL,
                   DebugLevel="Normal", 
                   iniFile=NULL )
{
# old default: 
#    Xmin=-0.5*(XSamples-1), 
#    Ymin=-0.5*(YSamples-1), 
#    DeltaX=1, 
#    DeltaY=1,
#    Xmin=-(XSamples)*sqrt(1/2), 
#    Ymin=-(YSamples)*sqrt(1/2), 
#    DeltaX=abs(2*(Xmin)/XSamples), 
#    DeltaY=abs(2*(Ymin)/YSamples),
#

# ========================================================================
#
# Description:
# NOTE: For each parameter the data-type is specified. In R you can use the
#       normal double-format for integer values, but note they become 
#       truncating to integer.
# NOTE: The C-routines using the float-format. The double precision-values
#       are  appropriate converted in C.
#
# This programm call a C-routine, which contains routines to call the
# different direct reconstruction routines. The mode of operation
# in determined by the 'mode'-
#
# iradon <- function(rData, XSamples, YSamples, mode="BF", 
#                    Xmin=-170, Ymin=-170, DeltaX=1.35, DeltaY=1.35,
#                    DebugLevel="Normal", oData=NULL )
#
# rData    (matrix) A matrix, that contain the sinogram image (for the 
#          reconstruction functions). 
# XSamples (integer) Number of samples on the x-axis in reconstructed image.
# YSamples (integer) Number of samples on the y-axis in reconstructed image.
# mode     (character) The function the program should do. Currently supported
#          functions are
#  "FB"    Filtered Backprojection.
#  "BF"    Filtering After Backprojection.
#  "CNF"   Central Slice. FFT based with Nearest Neighbor approximation.
#  "CBF"   Central Slice. FFT based with Bilinear Interpolation.
#  "CC"    Using a variation of the Central Slice. By the use of nonlinear 
#          sampling of the Radon domain with the use of the Chirp-z algorithm, 
#          the interpolation can be reduced to simple one dimensional linear 
#          interpolation compared to normal CS.
#  "CNC"   Central Slice. Chirp-z based with Nearest Neighbor approximation.
#  "CBC"   Central Slice. Chirp-z based with Bilinear Interpolation.
#  "Test"  Mode to run all the avaliable reconstruction routines on the same
#          image and compare with a L1Norm and L2Norm. Surrender value is "NULL".      
# Xmin     (double) The minimum x-position of the reconstructed image.
# Ymin     (double) The minimum y-position of the reconstructed image.
# DeltaX   (double) Sampling distance on the x-axis.
# DeltaY   (double) Sampling distance on the y-axis.
# InterPol    (integer) Interpolation level. Only used by Filtered Backprojection 
#             (FB).
# FilterTyp   (character) Filter type. Only used by Filtered Backprojection (FB).
#             Default: "Hamming1".
#  "Ramp"        Very good results by data without noise.
#  "Hamming1"    Generalized Hamming Filter with alpha=0.5 
#  "Hamming2"    Generalized Hamming Filter with alpha=0.54 
#                Both, Hamming1 and Hamming2 Better than "Ramp" in case of more 
#                noisy data.
# DebugLevel  (character) This parameter controls the level of output. The 
#             parameter is mixed and overrules with the one used in the
#             Print-statements.
#  "Normal"   Standard level of output to screen.
#  "Detail"   Allmost all output is logged to screen.
#  "HardCore" No information at all.
# oData    If mode='Test', measures of misfit can be made between the rData
#          and the oData. 'oData' contain a matrix with orginal data of the
#          radon tranformed data.
#
# ========================================================================

 ################################################################
 #       If an ini-file is specified then is reading now
 #
  args <- match.call()
  if (is.character(iniFile)){
      iniList <- readIni(iniFile, DebugLevel=DebugLevel)
  
      mode <- iniList$mode
      rData <- iniList$rData
      XSamples <- iniList$XSamples
      YSamples <- iniList$YSamples
      Xmin <- iniList$Xmin
      Ymin <- iniList$Ymin
      DeltaX <- iniList$DeltaX
      DeltaY <- iniList$DeltaY
      InterPol <- iniList$InterPol
      FilterTyp <- iniList$FilterTyp
      oData <- iniList$oData
      DebugLevel <- iniList$DebugLevel


  } else if (!is.null(iniFile))
      stop("'iniFile' must specified an INI-file or set to FALSE.")

 ################################################################
 #              Checking the input-parameter
 #
  if (!(is.matrix(rData)))
      stop("'rData' has to be of type 'matrix'.")

  if (as.integer(XSamples)!=XSamples){
      cat("WARNING : XSamples is not of type integer and is truncated to", (as.integer(XSamples)), "\n")
      }
  if (as.integer(YSamples)!=YSamples){
      cat("WARNING : YSamples is not of type integer and is truncated to", (as.integer(YSamples)), "\n")
      }
  if (as.integer(InterPol)!=InterPol){
      cat("WARNING : InterPol is not of type integer and is truncated to", (as.integer(InterPol)), "\n")
      }
  if (!(FilterTyp=="Ramp" || FilterTyp=="Hamming1" || FilterTyp=="Hamming2")){
      cat("WARNING: FilterTyp=''",FilterTyp,"'' is not supported. \n", sep="")
      cat("Default is used: FilterTyp=''Hamming1'' \n")
      FilterTyp<-"Hamming1"
      }
  if (!(DebugLevel=="Normal" || DebugLevel=="Detail" || DebugLevel=="HardCore")){
      cat("WARNING: DebugLevel=''",DebugLevel,"'' is not supported. \n", sep="")
      cat("Default is used: DebugLevel=''Normal'' \n")
      DebugLevel<-"Normal"
      }
    
  if (FilterTyp=="Ramp")
      FilterTypC <- "Ramp"
  else if (FilterTyp=="Hamming1")
      FilterTypC <- "Hanning"
  else if (FilterTyp=="Hamming2")
      FilterTypC <- "Hamming"

  if (mode=="CC" || mode=="CNC" || mode=="CBC" || mode=="CNF" || 
      mode=="CBF" || mode=="FB" || mode=="BF" ){

      rDataDim  <- dim(rData)
    
     # calling the C-routine "iradon" from iradon.c
      irData <-.C("iradon", 
                 as.double(rData), 
                 irData=double(XSamples*YSamples), 
                 as.character(mode),
                 as.integer(InterPol), 
                 as.character(FilterTypC),
                 as.character(DebugLevel), 
                 as.double(Xmin), 
                 as.double(Ymin), 
                 as.double(DeltaX), 
                 as.double(DeltaY), 
                 as.integer(rDataDim[1]), 
                 as.integer(rDataDim[2]), 
                 as.integer(XSamples), 
                 as.integer(YSamples),
                 PACKAGE="PET")$irData
    
      irData <- scaleImage(matrix(irData,nrow=XSamples,ncol=YSamples, byrow=TRUE))
      z <- list(irData=irData, 
            Header=list(SignalDim=c(XSamples,YSamples), 
                        XYmin=c(Xmin, Ymin), 
                        DeltaXY=c(DeltaX,DeltaY)), 
            call=args)
      class(z) <- "pet"

  } else if(mode=="Test"){
      if(!(is.matrix(oData)))
          stop("oData has to be of type 'matrix'.")    
      rDataDim  <- dim(rData)
      irDataDim <- dim(oData);
      if ( any(irDataDim != c(XSamples, YSamples)) ){
          cat("WARNING: If mode='Test' then the dimension of the inverse radon data has  to be the same as the dimension of the orginal data. \n")
          cat("Defaults are used: ")
          cat("XSamples=",irDataDim[1],", YSamples=",irDataDim[2],".\n",sep="")
      }
     # calling the C-routine "iradon" from iradon.c
      irData <-.C("iradon", 
                 as.double(rData), 
                 irData=as.double(oData), 
                 as.character(mode),
                 as.integer(InterPol), 
                 as.character(FilterTypC),   
                 as.character(DebugLevel), 
                 as.double(Xmin), 
                 as.double(Ymin), 
                 as.double(DeltaX), 
                 as.double(DeltaY), 
                 as.integer(rDataDim[1]), 
                 as.integer(rDataDim[2]), 
                 as.integer(irDataDim[1]), 
                 as.integer(irDataDim[2]),
                 PACKAGE="PET")$irData
      z <- NULL
  } else
      stop("The mode=",mode," is not supported.")

  return(z)
}
