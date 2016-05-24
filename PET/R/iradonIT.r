#########################################################################
#
# Copyright Joern Schulz,
#           Humboldt Universitaet zu Berlin, 
#           Germany 2006
# *********************************************************
#
# Name:          iradonIT.r
#                ---------------
# Author:        Joern Schulz
# Stand:         08.08.2006
#
# Uses:          
#
#########################################################################
# 
#
iradonIT <- function(rData, XSamples, YSamples,
                     StartImage="None",
                     mode="EM",
                     UseFast=1, 
                     RadonKernel="NN",
                     Iterations=20, 
                     IterationsType="random",
                     SaveIterations=0,
                     SaveIterationsName="",
                     LowestALevel=0.0,
                     ConstrainMin=-1,
                     ConstrainMax=-1,
                     Alpha=1.0, 
                     Beta=1.0,
                     Regularization=0,
                     KernelFileSave=0,
                     KernelFileName="",
                     RefFileName="None",
                     ThetaSamples=nrow(rData),
                     RhoSamples=ncol(rData),
                     ThetaMin=0,
                     RhoMin=-0.5*((2*round(sqrt(XSamples^2+YSamples^2)/2)+1)-1),
                     DeltaTheta=pi/ThetaSamples,
                     DeltaRho=(2*abs(RhoMin)+1)/RhoSamples,
                     Xmin=-0.5*(XSamples-1),
                     Ymin=-0.5*(YSamples-1),
                     DeltaX=1,
                     DeltaY=1,
                     OverSamp=0,
                     DebugLevel="Normal",
                     iniFile=NULL )
{
# old default: 
#    RhoMin=-0.5*(RhoSamples-1),
#    DeltaRho=1,
#    RhoMin=-1,
#    DeltaRho=2/RhoSamples,
#    Xmin=-1,
#    Ymin=-1,
#    DeltaX=2/XSamples,
#    DeltaY=2/YSamples,

# ========================================================================
#
# Description:
# This programm call a C-routine 'it.c', that contain routines to call
# different  iterative reconstruction routines. The mode of operation
# in determined by the 'mode'.
# NOTE: For each parameter the data-type is specified. In R you can use the
#       normal double-format for integer values, but note they become 
#       truncating to integer.
# NOTE: The C-routines using the float-format. The double precision-values
#       are  appropriate converted in C.
#
# iradonIT <- function(rData, XSamples, YSamples, mode="EM", UseFast=1, 
#     RadonKernel="NN", Iterations=20, IterationsType="random",
#     SaveIterations=0,  SaveIterationsName="", LowestALevel=0.0,
#     ConstrainMin=0, ConstrainMax=-1, Alpha=1.0, Beta=1.0,
#     Regularization=0, KernelFileSave=0, KernelFileName="",
#     ThetaSamples=nrow(rData), RhoSamples=ncol(rData), ThetaMin=0,
#     RhoMin=-1, DeltaTheta=pi/ThetaSamples, DeltaRho=2/(RhoSamples-1),
#     Xmin=-1, Ymin=-1, DeltaX=2/(XSamples-1), DeltaY=2/(YSamples-1),
#     DebugLevel="Normal" )
#
# rData    (matrix) A matrix, that contain the sinogram image (for the
#          reconstruction functions). The rows of the image correspond
#          to the sampled angles 'theta' and the columns correspond to 
#          the sampled distance. This both parameter determines the lines.
#          NOTE: For default, it is assumed that the sinogram is standardised
#          to a coordinate system with 'Theta' from 0 to pi and 'Rho' from
#          -1 to 1. Also it is assumed that the reconstructed image is
#          standardised to a coordinate system with 'XSamples' from -1 to 1 
#          and 'YSamples' from -1 to 1.
# XSamples (integer) Number of samples on the x-axis in reconstructed image.
# YSamples (integer) Number of samples on the y-axis in reconstructed image.
# StartImage   (matrix) If provided the initial guess on the reconstructed 
#          image is taken from this image, else a constant solution will be
#          assumed from start. Default: StartImage="None"
# mode     (character) The iterative reconstruction method that the program
#          should do.
#          Default: mode="EM". Currently supported functions are
#  "ART"   Algebraic Reconstruction Technique
#  "EM"    Likelihood Reconstruction using Expectation Maximation
#  "CG"    Least Squares Conjugate Method
# UseFast  (integer) If 1 then fast reconstruction is used, where the system
#          matrix is stored using sparse techniques, else slower but memory
#          efficient reconstruction is used. Default: UseFast=1.
# RadonKernel   (character) The type of kernel used to model the system matrix. 
#          Default: RadonKernel="NN". Currently available are the following    
#          methods:
#  "NN"    Two-level Nearest Neighbor approximation. (Memory consuming).
#  "RNN"   Ray driven Nearest Neighbor discrete Radon transform based (Very 
#          fast with small system matrix).
#  "RL"    Ray driven Linear Interpolation discrete Radon transform based 
#          (Fast with small system matrix).
#  "P1"    Method based on Radon transformation of square with pre-guidance 
#          (slow but good).
#  "P2"    Method based on Radon transformation of square with no pre-guidance 
#          (slower but better).
#  "SINC"  Sinc interpolation methods in the image and analytically Radon of
#          that. (Very slow). 
# Iterations   (integer) For EM and CG the number of iterations before the
#          iteration ends.
#          For ART the number of full iterations, i.e., divided by the number 
#          of rows in the system matrix.  Default: Iterations=20. 
# IterationsType   (character) Two iterations type for ART are possible:
#          "cyclic" selection of the row index or "random" selection. 
#          Default: IterationsType="random".
# SaveIterations   (integer) If set to 1 the current solution will be saved
#          after each iteration under the name
#          SaveIterationsName+CurrentIteration (also possible 2,3 ...).
# SaveIterationsName    (character) In Case of 'SaveIteration' is non-zero., 
#          the parameter determined the name of currently saved iterations.
#          It contains the path to filename and the filename. The path
#          containing the path relatively to the working-directory of your
#          R-session or containig the full path to the file.
#          Default: SaveIterationsName="".
# LowestALevel   (double) If fast reconstruction is used, the the matrix
#          elements are truncated to this level relative to the sampling 
#          distance of x. Default: LowestALevel=0.0.
# ConstrainMin   (double) After each iteration the solution in each sample
#          will forced above this limit.  Default: ConstrainMin=0.
# ConstrainMax   (double) After each iteration the solution in each sample
#          will forced below this limit. Default: ConstrainMax=-1.
#          NOTE: For both limit it is assumed that negative limits imply that 
#          the feature is not used.
# Alpha    (double) For ART the initial update weight (if not specified then it
#          is set to 1). Default: Alpha=1.0.
# Beta     (double) For ART the multiplicative change to the weight factor,
#          which should be less than one (if not specified then it is set to 1). 
#          Default: Beta=1.0.
# Regularization   (double) If set to 1 and using fast reconstruction, rows will
#          be appended to the system matrix with a simple Laplace operator. A
#          weight factor should also be incorporated. 
#          Default: Regularization=0.
# KernelFileSave    (integer) If 1 then the system matrix is saved under the
#          the name 'KernelFileName'. Default: KernelFileSave=0.
# KernelFileName   (character) If using fast reconstruction, the system matrix
#          will be saved and restored with this sif-name. If the system matrix
#          read is incopatible with the sampling parameters, then a new will be
#          generated. KernelFileName contains the path to filename and the
#          filename. The path containing the path relatively to the
#          working-directory of your R-session or containig the full path 
#          to the file. Default: KernelFileName="".
# RefFileName   (character) Error measures can be made between the image 
#          containing in the RefFileName and reconstructed image from \code{rData}. 
#          Defaults to \code{RefFileName="None"}.
# ThetaSamples   (integer) Number of angular samples $T$ in the sinogram. 
#          Default: ThetaSamples=nrow(rData).
# RhoSamples   (integer) Number of samples in the sinogram $R$ in the
#          $\rho$-direction. Default: RhoSamples=ncol(rData).
# ThetaMin (double) Start of the angular sampling (should be set to zero.). 
#          Default: ThetaMin=0.
# RhoMin   (double) Start of sampling positions in $\rho$ (should be set to
#          $-\Delta\rho\frac{R-1}{2}$.  Default: RhoMin=-1.
# DeltaTheta   (double) Angular sampling distance (should be set to
#          pi/ThetaSamples). Default: DeltaTheta=pi/ThetaSamples.
# DeltaRho (double) Sampling distance in $\rho$, i.e., $\Delta \rho$. 
#          Default: DeltaRho=2/(RhoSamples-1).
# Xmin     (double) The minimum x-position of the reconstructed image.
#          Default: Xmin=-1.
# Ymin     (double) The minimum y-position of the reconstructed image.
#          Default: Ymin=-1.
# DeltaX   (double) Sampling distance on the x-axis.
#          Default: DeltaX=2/(XSamples-1).
# DeltaY   (double) Sampling distance on the y-axis.
#          Default: DeltaY=2/(YSamples-1).
# OverSamp   (integer)  Use oversampling - Using squared number of samples. 
#          Defaults to \code{OverSamp=0}. }
# DebugLevel   (character) This parameter controls the level of output. 
#          Default: DebugLevel="Normal".
#          Three possibilties are available:
#   "Detail"   Detail information output to screnn.
#   "Normal"   Standard level of output to screen.
#   "HardCore" No information.
# iniFile   (character)  If \code{iniFile!=NULL}, then \code{iniFile} have to 
#          be the name of an ini-file, including phatname to the file. In the 
#          case of specified \code{iniFile} all parameter are reading from the file.
#
#
#########################################################################

  args <- match.call()
  MaxSaves <- 200 # limit for Iterations/SaveIterations > MaxSaves

 ################################################################
 #       If an ini-file is specified then is reading now
 #
  if (is.character(iniFile)){
      iniList <- readIni(iniFile, DebugLevel=DebugLevel)
    
      mode <- iniList$mode
      rData <- iniList$rData
      XSamples <- iniList$XSamples
      YSamples <- iniList$YSamples
      StartImage <- iniList$StartImage
      UseFast <- iniList$UseFast
      RadonKernel <- iniList$RadonKernel
      Iterations <- iniList$Iterations
      IterationsType <- iniList$IterationsType
      SaveIterations <- iniList$SaveIterations
      SaveIterationsName <- iniList$SaveIterationsName
      LowestALevel <- iniList$LowestALevel
      ConstrainMin <- iniList$ConstrainMin
      ConstrainMax <- iniList$ConstrainMax
      Alpha <- iniList$Alpha
      Beta <- iniList$Beta
      Regularization <- iniList$Regularization
      KernelFileSave <- iniList$KernelFileSave
      KernelFileName <- iniList$KernelFileName
      RefFileName <- iniList$RefFileName
      ThetaSamples <- iniList$ThetaSamples
      RhoSamples <- iniList$RhoSamples
      ThetaMin <- iniList$ThetaMin
      RhoMin <- iniList$RhoMin
      DeltaTheta <- iniList$DeltaTheta
      DeltaRho <- iniList$DeltaRho
      Xmin <- iniList$Xmin
      Ymin <- iniList$Ymin
      DeltaX <- iniList$DeltaX
      DeltaY <- iniList$DeltaY
      OverSamp <- iniList$OverSamp
      DebugLevel <- iniList$DebugLevel

  } else if (!is.null(iniFile))
      stop("'iniFile' must specified an INI-file or set to FALSE.")

 ################################################################
 #              Checking the input-parameter
 #
  if (!(is.matrix(rData)))
      stop("'rData' has to be of type 'matrix'.")

  if (!is.character(StartImage)){
      if (!(is.matrix(StartImage)))
          stop("'StartImage' has to be of type 'matrix'.")
      if (nrow(StartImage)!=XSamples){
          cat("WARNING: 'StartImage' has to be the same numbers of rows as the reconstructed image, specified by 'XSamples'. \n")
          cat("Default is used: StartImage='None'\n")
          StartImage<-"None"
      }
      if (ncol(StartImage)!=YSamples){
          cat("WARNING: 'StartImage' has to be the same numbers of columns as the reconstructed image, specified by 'YSamples'. \n")
          cat("Default is used: StartImage='None'\n")
          StartImage<-"None"
      }
  } else if (StartImage!="None"){
      cat("WARNING: 'StartImage' has to be 'None' or of type 'matrix'.\n")
      cat("StartImage='", StartImage,"'  is not supported. Default is used. \n", sep="")
  }

  if (!is.character(StartImage)) StartImageTrue <- 1
  else {                  StartImageTrue <- 0
                          StartImage <- 0     }

  if (as.integer(XSamples)!=XSamples){
      cat("WARNING : 'XSamples' is not of type integer and is truncated to", (as.integer(XSamples)), "\n")
  }
  if (as.integer(YSamples)!=YSamples){
      cat("WARNING : 'YSamples' is not of type integer and is truncated to", (as.integer(YSamples)), "\n")
  }
  if (as.integer(Iterations)!=Iterations){
      cat("WARNING : 'Iterations' is not of type integer and is truncated to", (as.integer(Iterations)), "\n")
  }
  if (as.integer(ThetaSamples)!=ThetaSamples){
      cat("WARNING : 'ThetaSamples' is not of type integer and is truncated to", (as.integer(ThetaSamples)), "\n")
  }
  if (as.integer(RhoSamples)!=RhoSamples){
      cat("WARNING : 'RhoSamples' is not of type integer and is truncated to", (as.integer(RhoSamples)), "\n")
  }
  if (as.integer(OverSamp)!=OverSamp){
      cat("WARNING : 'OverSamp' is not of type integer and is truncated to", (as.integer(OverSamp)), "\n")
  }
    
  if (!(mode=="EM" || mode=="ART" || mode=="CG")){
      stop("Unknown iterative reconstruction method: 'mode'=",mode)
  }
  if (!(RadonKernel %in% c("NN", "RNN", "RL", "P1", "P2"))){
      cat("WARNING: Kernel of type '",RadonKernel,"' is not supported. \n",sep="")
      cat("Default is used: RadonKernel='NN' \n")
      RadonKernel<-"NN"
  }
  if (!(IterationsType=="random" || IterationsType=="cyclic")){
      cat("WARNING: IterationsType='",IterationsType,"' is not supported. \n", sep="")
      cat("Default is used: IterationsType='random' \n")
      IterationsType<-"random"
  }
  if (!(DebugLevel=="Normal" || DebugLevel=="Detail" || DebugLevel=="HardCore")){
      cat("WARNING: DebugLevel='",DebugLevel,"' is not supported. \n", sep="")
      cat("Default is used: DebugLevel='Normal' \n")
      DebugLevel<-"Normal"
  }
    
  if (UseFast!=0 && UseFast!=1){
      cat("WARNING: UseFast should be only set to 0 for FALSE or 1 for TRUE. \n")
      cat("Default is used: UseFast=1 \n")
      UseFast<-1
  }
  if (SaveIterations && Iterations/SaveIterations > MaxSaves){
      cat("WARNING: To large number of saves requested, aborting SaveIterations. \n")
      cat("Default is used: SaveIterations=0 \n")
      SaveIterations<-0
  }
  if (KernelFileSave!=0 && KernelFileSave!=1){
      cat("WARNING: KernelFileSave should be only set to 0 for FALSE or 1 for TRUE. \n")
      cat("Default is used: KernelFileSave=0 \n")
      KernelFileSave<-0
  }

  if (!(is.character(SaveIterationsName))){
      stop("SaveIterationsName have to be of type character. It contains the path to filename and the filename. The path containing the path relatively to the working-directory of your R-session or containig the full path to the file.")
  }  
  if (SaveIterationsName!=""){
      DirName <- dirname(SaveIterationsName)
      if (file.access(DirName,0)){
          cat("The directory '", DirName,"' , where the iterations should be saving, doesn't exist. \n", sep="")
          stop()
  }}

  if (!(is.character(KernelFileName))){
      stop("KernelFileName have to be of type character. It contains the path to filename and the filename. The path containing the path relatively to the working-directory of your R-session or containig the full path to the file.")
  }  
  if (KernelFileName!=""){
      DirName <- dirname(KernelFileName)
      if (file.access(DirName,0)){
          cat("The directory '", DirName,"' for saving the system matrix doesn't exist. \n", sep="")
          stop()
  }}

  if (!(is.character(RefFileName))){
      stop("RefFileName have to be of type character. It contains the path to filename and the filename. The path containing the path relatively to the working-directory of your R-session or containig the full path to the file.")

  } else if (RefFileName!="None"){
      if (file.access(RefFileName,0)){
          cat("The adress '", RefFileName,"' for RefFileName doesn't exist. \n", sep="")
          stop()
      } 
  } else RefFileName <- ""

  if (SaveIterationsName=="" && SaveIterations==1){
      cat("WARNING: It is not possible to save iterations, because SaveIterationsName='' and therefore it is not specified. \n")
      cat("Default is used: SaveIterations=0 \n")
      SaveIterations<-0
  }
  if (KernelFileName=="" && KernelFileSave==1){
      cat("WARNING: It is not possible to save the system matrix, because         KernelFileName='' and therefore it is not specified.  \n")
      cat("Default is used: KernelFileSave=0 \n")
      KernelFileSave<-0
  }
    
 
 # --------->>>>>>>>>>> current warnings <<<<<<<<<<<<<<<-----------
  if(!(DebugLevel=="HardCore")){
  if (UseFast==0){
      if (RadonKernel=="RNN" || RadonKernel=="RL"){
          cat("\n")
          cat("WARNING: At the moment, in case of a slow iterative reconstruction method and RadonKernel=='RNN' or 'RL' the computation of the system matrix is probably not correct! The mistake is probably in the C-Code of 'amatrix.c' in the procedure 'GenerateAMatrixColumn'. I invite you to find the mistake and send me an Email ;o). \n")
          cat("\n")
      }
  }
  if (mode=="CG"){
      cat("\n")
      cat("WARNING: At the moment, it cannot be ensured that the CG algorithm work correct. I invite you to look into the C-code of 'cg.c'. The problem is described there in the procedures 'FAST_CG' and  'SLOW_CG'. \n") 
      cat("\n")
  }
  }
 
 # ----------------------------------------------------------------------
 # calling the C-routine "it" from it.c
 #

  ITradonData <- .C("it", 
            as.double(rData), 
            ITradonData=double(XSamples*YSamples),
            as.integer(StartImageTrue),
            as.double(StartImage),
            as.character(mode),
            as.integer(UseFast),
            as.character(RadonKernel), 
            as.character(IterationsType),
            as.integer(Iterations),
            as.integer(SaveIterations),
            as.character(SaveIterationsName),
            as.double(LowestALevel), 
            as.double(ConstrainMin), 
            as.double(ConstrainMax), 
            as.double(Alpha),
            as.double(Beta),
            as.double(Regularization),
            as.integer(KernelFileSave),
            as.character(KernelFileName),
            as.character(RefFileName),
            as.integer(ThetaSamples),
            as.double(ThetaMin),
            as.double(DeltaTheta),
            as.integer(RhoSamples),
            as.double(RhoMin),
            as.double(DeltaRho),
            as.double(Xmin),
            as.double(Ymin),
            as.double(DeltaX),
            as.double(DeltaY),
            as.integer(XSamples),
            as.integer(YSamples),
            as.integer(OverSamp),
            as.character(DebugLevel),
            PACKAGE="PET")$ITradonData
    
  ITradonData <- scaleImage(matrix(ITradonData,nrow=XSamples,ncol=YSamples, byrow=TRUE))

  z <- list(irData=ITradonData, 
            Header=list(SignalDim=c(XSamples,YSamples), 
                       XYmin=c(Xmin, Ymin), 
                       DeltaXY=c(DeltaX,DeltaY)), 
            call=args)
  class(z) <- "pet"
  return(z)
}
