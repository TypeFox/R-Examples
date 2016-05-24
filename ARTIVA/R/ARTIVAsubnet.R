ARTIVAsubnet <-
function(targetData, parentData, targetName="Target",parentNames=NULL, 
           dataDescription=NULL, saveEstimations=TRUE, saveIterations=FALSE, 
           savePictures = TRUE, outputPath=NULL, dyn=1,  segMinLength=2, maxCP=NULL, 
           maxPred=NULL, nbCPinit=NULL, CPinit=NULL, niter=50000, burn_in=NULL, 
           PSRFactor=FALSE, PSRF_thres=1.05, segmentAnalysis=TRUE,  edgesThreshold=0.5,
           layout = "fruchterman.reingold", cCP= 0.5, cEdges=0.5, alphaCP=1,
           betaCP=0.5, alphaEdges=1, betaEdges=0.5, v0=1, gamma0=0.1, alphad2=2,
           betad2=0.2, silent=FALSE){
 pdf.options(fonts="serif")

  # data must not contain missing values
  if(sum(is.na(targetData))>0){
   stop("\n*****************************\n
ERROR: targetData must not contain missing values \n
*****************************\n")
 }
   if(sum(is.na(parentData))>0){
   stop("\n*****************************\n
ERROR: parentData must not contain missing values \n
*****************************\n")
  }
  
  # targetData must be a vector
  if(!is.vector(targetData)){
    stop("\n*****************************\n
ERROR: targetData must be a numeric vector\n
*****************************\n")
  }

  # parentData must be a matrix or a vector
  if(!is.matrix(parentData)){
    if(!is.vector(parentData)){
      stop("\n*****************************\n
ERROR: parentData must be a matrix (if several parent genes are proposed) or a numeric vector (if only one parent gene is proposed) \n
*****************************\n")
      }
    }  

  # If parentData is a vector, it is converted into a matrix
  if(!is.matrix(parentData)){
    parentData=t(as.matrix(parentData))
  }

  # targetData and parentData must have the same number of expression measurements
  if(length(targetData) != ncol(parentData)){
    stop("\n*****************************\n
ERROR: targetData and parentData don't have the same number 
of expression measurements\n
*****************************\n")
  }
 
  # test of the vector describing repeated experiments
  if(!(is.null(dataDescription))){
    if(length(unique(table(dataDescription)))>1){
      stop("\n*****************************\n
ERROR: DataDescription is not correctly writen. 
Please give the same number of repetitions for each time point measurement\n
*****************************\n")
    }
  }
  # test savePictures & saveEstimations
  if(!silent & !savePictures){
    print("--- !!! WARNING MESSAGE !!! ----")
    print("Plots will not be saved (see the parameter savePictures)")
    print("---------------------------------")
  }
  if(!silent & !saveEstimations){
    print("--- !!! WARNING MESSAGE !!! ----")
    print("Estimation values will not be saved (see the parameter saveEstimations)")
    print("---------------------------------")
  }
  
  # number of timepoints (n) and repetitions (m)
  if(is.null(dataDescription)){
    n = length(targetData)
    
    m = 1
  }else{
    n=length(unique(dataDescription))
    m=length(targetData)/n
  }

  # Position of each time point in the data (designed for the algorithm)
  Mphase=seq(1,n*m+1,by=m)-dyn*m

  # number of parent genes
  q=dim(parentData)[1]

  # maximal number of changepoints (CPs)
  if(is.null(maxCP)){
    maxCP=round(min((n-1-dyn)/segMinLength-1,15))
  }

  #  number of changepoints (CPs) at initialization
  if(is.null(nbCPinit)){
    nbCPinit=sample(0:round((maxCP/2)),1)  #min(floor(n/2),3, maxCP)
  }

  # maximal number of Factors (TFs)	
  if(is.null(maxPred)){
	maxPred=min(dim(parentData)[1],15)
  }
  
  # number of iterations used for the output=>> postDist functions.
  if(is.null(burn_in))burn_in=round(niter/4)

  # nbVarMax = maximal number of variances (default: kmax, put it to 1 if the variance is the same for all phases)
  nbVarMax=maxCP+1
     
  # names of predictors
  if(is.null(parentNames)){
    parentNames=as.character(paste("pa_", c(1:nrow(parentData)), sep = ""))
  }

  # values for PSRF factor computation (if wanted)
  # PSRF computation frequency (number of iterations) after the minimal burn-in period 
  PSRFactor_freq=min(niter,1000)
  # Number of sequences in which the second half of the iterations is devided for PSRF computation
  PSRFactor_nbSeq=10
  # Number of extra iterations after convergence is established (these samples will be used for estimating the posterior distribution)
  PSRFactor_nbIterEval=min(niter,50000)

  ## General information concerning the analyzed data set
  if(!silent)
  {
    print("========================")
    print("GENERAL INFORMATION")
    print(paste(" *** Name of the analyzed target gene:", targetName))
    print(paste(" *** Number of different time point measurements:", n))
    print(paste(" *** Number of repeated measurements:", m))
    print(paste(" *** Number of potential parent gene(s):", q))
    print("========================")
    print("ARTIVA PARAMETERS")
    print(paste(" **** Time delay considered in the auto-regressive process:", dyn))
    print(paste(" **** Minimal length to define a temporal segment:", segMinLength, "time points"))
    print(paste(" **** Maximal number of changepoints (CPs):", maxCP))
    print(paste(" **** Number of CPs at the algorithm initialization:", nbCPinit))
    if(!is.null(CPinit)){
      print(paste(" **** Initial positions of CPs at the algorithm initialization:", nbCPinit))
    }
    print(paste(" **** Number of iterations:", niter))
    print(paste(" **** Is PSRF factor calculated?", PSRFactor))
    print("========================")
  }

  # Create Global Variables used in all functions
  GLOBvar = list(n=n, p=1, q=q, qmax=maxPred, smax=maxCP, dyn=dyn, segMinLength=segMinLength, nbVarMax=nbVarMax, m=m, Mphase=Mphase, niter=niter, targetName=targetName, parentNames=parentNames, silent=silent )

  # Create HyperParms Variables used in all functions
  HYPERvar = list(cD=cCP, alphaD=alphaCP, betaD=betaCP, c=cEdges, v0=v0, gamma0=gamma0, alphad2=alphad2, betad2=betad2, alphalbd=alphaEdges, betalbd=betaEdges)

  # Create Output Variables used in output functions
  OUTvar = list(burn_in=burn_in, PSRFactor=PSRFactor, PSRFactor_freq=PSRFactor_freq,PSRFactor_nbSeq=PSRFactor_nbSeq, PSRFactor_nbIterEval=PSRFactor_nbIterEval,PSRF_thres=PSRF_thres)
  
 
  # build response Y and predictor X
  input = buildXY(targetData, parentData, dataDescription, m, dyn)
  X = input$X
  Y = input$Y

  # initialize system
  if(!silent)print("STEP 1: Starting the initialisation procedure")
  initiation = init(X, nbCPinit, GLOBvar, HYPERvar, OUTvar, CPinit)

  # run niter iterations
  if(!silent)print(paste("STEP 2: Starting the ARTIVA RJ-MCMC procedure"))
  runiteration = main(X, Y, initiation, GLOBvar, HYPERvar, OUTvar)

  initPath=getwd()
  if(saveEstimations | saveIterations | savePictures){
    
    if(is.null(outputPath)){
      # Output directory path (if it doesn t exist create it)
      if(.Platform$OS.type == "unix"){
        if(! "ARTIVAsubnet" %in% system("ls" ,intern=TRUE)) {
          system("mkdir ARTIVAsubnet", ignore.stderr = TRUE)
        }  
      }else{# if(.Platform$OS.type == "unix"){
        shell("mkdir ARTIVAsubnet", intern = TRUE,mustWork =NA)
      }
      outputPath="ARTIVAsubnet"
      
    }else{#if(is.null(outputPath)){
      testPath=strsplit(outputPath, split="/")[[1]]
      if(length(testPath)>1){
        setwd(substr(outputPath,0,nchar(outputPath)-nchar(testPath[length(testPath)])-1))
        outputPath=strsplit(outputPath, split="/")[[1]][length(strsplit(outputPath, split="/")[[1]])]    
      }
      
      if(.Platform$OS.type == "unix"){
        if(! outputPath  %in% system( "ls ",intern=TRUE)) {
          system(paste("mkdir ", outputPath, sep=""))
        }
      }else{# if(.Platform$OS.type == "unix"){
        shell(paste("mkdir ", outputPath, sep=""), intern = TRUE,mustWork =NA)
      }
    }
      ## Output Stocks directory path (to store the other results of the procedure) (if it doesn t exist create it)
      if(saveIterations){
        if(.Platform$OS.type == "unix"){
         if(! "Iterations" %in% system(paste( "ls ",outputPath),intern=TRUE)) {
           system(paste("mkdir ",outputPath,"/Iterations",sep=""))
         }
       }else{# if(.Platform$OS.type == "unix"){
         shell(paste("mkdir ",outputPath,"/Iterations",sep=""), intern = TRUE,mustWork =NA)
       }
       
      outputStockPath=paste(outputPath,"/Iterations",sep="")

      if(!silent)print(paste("Saving Iterations at: ", getwd(), "/",outputStockPath),sep="")
     
      write.table(as.matrix(runiteration$samples$CP),paste(outputStockPath,"CPsamples_target_",targetName,".txt",sep=""))
      write.table(as.matrix(runiteration$samples$Edges),paste(outputStockPath,"EdgesSamples_target_",targetName,".txt",sep=""))
      write.table(as.matrix(runiteration$samples$coeff),paste(outputStockPath,"CoeffSamples_target_",targetName,".txt",sep=""))
      write.table(as.matrix(runiteration$samples$variance),paste(outputStockPath,"VarianceSamples_target_",targetName,".txt",sep=""))
    }
  }
    
  # CP analysis: compute the CP number and CP position posterior distribution 
   if(!silent)print("STEP 3: Computing posterior distributions for the changepoints (CPs) number and location")
  
  CPpostDist=CP.postDist(runiteration$samples$CP, burn_in=OUTvar$burn_in, segMinLength=GLOBvar$segMinLength)
 
  if(saveEstimations){
    if(.Platform$OS.type == "unix"){
      if(! "Estimations" %in% system(paste( "ls ",outputPath),intern=TRUE)) { system(paste("mkdir ",outputPath,"/Estimations",sep=""))}
    }else{# if(.Platform$OS.type == "unix"){
      tmp=getwd()
      setwd(outputPath)
      shell("mkdir Estimations", intern = TRUE,mustWork =NA)
      setwd(tmp)
    }
    outputResPath=paste(outputPath,"/Estimations",sep="")
    if(!silent)print(paste("Saving Estimations at: ", getwd(), "/",outputResPath,sep=""))
  
    write.table(as.matrix(CPpostDist$CPnumber),paste(outputResPath,"/CPnumberPostDist_",targetName,".txt",sep=""))
    write.table(as.matrix(CPpostDist$CPposition),paste(outputResPath,"/CPpositionPostDist_",targetName,".txt",sep=""))
    write.table(as.matrix(CPpostDist$CPpos),paste(outputResPath,"/estimatedCPposition_",targetName,".txt",sep=""))
  }

  segAnalysis=NULL
  if(segmentAnalysis){
    segAnalysis=ARTIVAsubnetAnalysis(targetData=targetData, parentData=parentData, CPpostDist=CPpostDist,CPsamples=runiteration$samples$CP,coefSamples=runiteration$samples$coef,TFnumber=GLOBvar$q,segMinLength=GLOBvar$segMinLength,edgesThreshold=edgesThreshold, burn_in=burn_in, targetName=targetName, parentNames=parentNames, CPpos=CPpostDist$estimatedCPpos,savePictures=savePictures, saveEstimations=saveEstimations,outputPath=outputPath, layout=layout, silent=silent, inARTIVAsubnet=TRUE, onepage=FALSE)
  }else{
    if(!silent)print(paste("SKIPPING STEP 4: the posterior distribution for the regulatory models in each temporal phase is NOT computed "))

  }

  
  if(!silent)print("...............Ending ARTIVA analysis.................")
  setwd(initPath)
  
  ## return all samples	
  return(list(Samples=runiteration$samples,Counters=runiteration$counters, CPpostDist=CPpostDist,  nbSegs=segAnalysis$nbSegs, SegmentPostDist=segAnalysis$SegmentPostDist, network = segAnalysis$network, GLOBvar=GLOBvar, HYPERvar=HYPERvar, OUTvar=OUTvar,targetData=targetData, parentData=parentData ))
}
