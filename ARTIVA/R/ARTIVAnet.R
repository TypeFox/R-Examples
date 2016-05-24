ARTIVAnet <-
function(targetData, parentData, targetNames = NULL, parentNames = NULL,
           dataDescription=NULL, saveEstimations=TRUE, saveIterations=FALSE,
           savePictures = TRUE, 
           outputPath=NULL, dyn=1,  segMinLength=2, 
           maxCP=NULL, maxPred=NULL, nbCPinit=NULL, 
           CPinit=NULL, niter=50000, burn_in=NULL, 
           PSRFactor=FALSE,PSRF_thres=1.05, 
           segmentAnalysis=TRUE, edgesThreshold=0.5, 
           layout = "fruchterman.reingold", cCP= 0.5, cEdges=0.5,
           alphaCP=1, betaCP=0.5, alphaEdges=1, betaEdges=0.5,
           v0=1, gamma0=0.1, alphad2=2, betad2=0.2, silent=FALSE)
{

  # targetData is a matrix with all the target gene that are successively 
  # analyzed with ARTIVA

  if(!silent){
    print("====================================================================")
    print("Starting  ARTIVAnet  procedure")
    print("====================================================================")
  }
  
  initPathARTIVAnet=getwd()

  if(is.null(outputPath)){
    ## Output directory path (if it doesn t exist create it)

	if(.Platform$OS.type == "unix"){
        if(! "ARTIVAnet" %in% system("ls" ,intern=TRUE))
		{
			system("mkdir ARTIVAnet")
		}
    }else{# if(.Platform$OS.type == "unix"){
      shell("mkdir ARTIVAnet", intern = TRUE, mustWork =NA)
    }	
	outputPath="ARTIVAnet"
	
  }else{
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
        shell(paste("mkdir ",outputPath, sep=""), intern = TRUE, mustWork =NA)
      }
	  
  }

  
  if(sum(is.na(parentData))>0){
    stop("\n**********************************************************\n
ERROR: parentData must not contain missing values \n
**********************************************************\n")
  }
  
  # targetData must be a matrix or a vector
  if(!is.matrix(targetData)){
    if(!is.vector(targetData)){
      stop("\n**********************************************************\n
ERROR:  targetData must be a matrix (if several target genes are proposed) or a numeric vector (if only one target gene is proposed) \n
**********************************************************\n")
      }
    }  

  # If parentData is a vector, it is converted into a matrix
  if(!is.matrix(targetData)){
    targetData=t(as.matrix(targetData))
  }

  # parentData must be a matrix or a vector
  if(!is.matrix(parentData)){
    if(!is.vector(parentData)){
      stop("\n**********************************************************\n
ERROR: parentData must be a matrix (if several parent genes are proposed) or a numeric vector (if only one parent gene is proposed) \n
**********************************************************\n")
      }
    }  

  # targetData and parentData must have the same number of expression measurements
 if(nrow(targetData)==1 || is.null(nrow(targetData))){
    nTargetTemp=length(targetData)
  }else{nTargetTemp=ncol(targetData)}

  if(nrow(parentData)==1 || is.null(nrow(parentData))){
    nParentTemp=length(parentData)
  }else{nParentTemp=ncol(parentData)}
  
  if(nTargetTemp != nParentTemp){
    stop("\n*****************************\n
ERROR: targetData and parentData don't have the same number 
of expression measurements\n
*****************************\n")
  }

  ## names of predictors
  if(is.null(targetNames)){
    targetNames=as.character(paste("tg_", c(1:nrow(targetData)), sep = ""))
  }
 
  
  # test of the vector describing repeated experiments
  if(!(is.null(dataDescription))){
    if(length(unique(table(dataDescription)))>1){
      stop("\n**********************************************************\n
ERROR: DataDescription is not correctly writen. 
Please give the same number of repetitions for each time point measurement\n
**********************************************************\n")
    }
  }
  
  # test savePictures & saveEstimations
  if(!savePictures){
    print("----------------------------------------------------------")
    print("---------------- !!! WARNING MESSAGE !!! -----------------")
    print("Plots will not be saved (see the parameter savePictures)")
    print("----------------------------------------------------------")
  }
  if(!saveEstimations){
    print("----------------------------------------------------------")
    print("---------------- !!! WARNING MESSAGE !!! -----------------")
    print("Estimation values will not be saved (see the parameter saveEstimations)")
    print("----------------------------------------------------------")
   }

  ## number of timepoints
  if(is.null(dataDescription)){
    n = dim(targetData)[2]
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
    maxCP=min((n-1-dyn)/segMinLength,15)
  }

  #  number of changepoints (CPs) at initialization
  if(is.null(nbCPinit)){
    nbCPinit= sample(0:round((maxCP/2)),1) #min(floor(n/2),5)
  }

  # maximal number of changepoints (CPs)
  if(is.null(maxCP)){
	maxCP=min(n-1-dyn,15)
  }

## General information concerning the analyzed data set
   if(!silent){
     print("====================================================================")
     print("GENERAL INFORMATION")
     print(paste(" *** Number of different time point measurements:", n))
     print(paste(" *** Number of repeated measurements:", m))
     print(paste(" *** Number of potential parent gene(s):", q))
     print("====================================================================")
     print("ARTIVA PARAMETERS")
     print(paste(" **** Time delay considered in the auto-regressive process:", dyn))
     print(paste(" **** Minimal length to define a temporal segment:", segMinLength, "time points"))
     print(paste(" **** Maximal number of changepoints (CPs):", maxCP))
     print(paste(" **** Number of CPs at the algorithm initialization:", nbCPinit))
     if(!is.null(CPinit)){
       print(paste(" **** Initial positions of CPs at the algorithm initialization:", CPinit))
     }
     print(paste(" **** Number of iterations:", niter))
     print(paste(" **** Is PSRF factor calculated?", PSRFactor))
     print("====================================================================")
    }

 
  GlobalNetwork = NULL
  

  estimationForTargets=NULL
  
  for(i in 1:nrow(targetData))
  {
     if(!silent){
       print(paste(" *** Name of the analyzed target gene:", targetNames[i]))
     }

     currentTargetData = targetData[i,]

    ## catch if error
    tryCatch({
      ARTIVAanalysis = ARTIVAsubnet(targetData = currentTargetData, 
				targetName = targetNames[i],
				parentData = parentData, 
				parentNames, 
				dataDescription, saveEstimations, 
				saveIterations, savePictures, 
				outputPath = outputPath, 
				dyn,  segMinLength, 
				maxCP, maxPred, nbCPinit, 
				CPinit, niter, burn_in, 
				PSRFactor,PSRF_thres, 
				segmentAnalysis, edgesThreshold, layout,
				cCP, cEdges, alphaCP, 
				betaCP, alphaEdges, betaEdges,
				v0, gamma0, alphad2, betad2, silent=TRUE)

      GlobalNetwork = rbind(GlobalNetwork, ARTIVAanalysis$network)
      estimationForTargets=c(estimationForTargets,i)
    }, error = function(ex) {
      print(ex);
    }) ## end catch if error

  }

  is.na(GlobalNetwork$PostProb) = 0
  
  ##save GlobalNetwork
  if(saveEstimations){
    write.table(GlobalNetwork, 
                  file = paste(outputPath,"/ARTIVA_FinalNetwork.txt",sep=""), 
                  sep = "\t", quote = FALSE)
  }
  
  if(savePictures)
    {
      pdf(paste(outputPath,"/ARTIVA_FinalNetwork.pdf",sep=""), width = 10, height = 7)
    }

  print("Calculation of graphs representations...")
  
  traceGeneProfiles(targetData, parentData)
  ##,CPpos = c(GlobalNetwork$CPstart, GlobalNetwork$CPend))
  if(savePictures)par(mfrow=c(1,1))
  traceNetworks(GlobalNetwork, edgesThreshold, layout=layout, onepage=!savePictures)

  #par(mfrow=c(1,1))
  geneNetworkSummary(GlobalNetwork, edgesThreshold)

  if(savePictures)
    {
      dev.off()
    }
   setwd(initPathARTIVAnet)
  
  estimationErrors = which( !(1:nrow(targetData) %in% estimationForTargets))
  if(length(estimationErrors)>0)
  {
    print("====================================================================")
    print( paste("WARNING :  No estimation for target gene:", estimationErrors))
    print("Please check the error.")
    print("====================================================================")
  }
  testPath=strsplit(outputPath, split="/")[[1]]
  if(length(testPath)==1){
    outputPath=paste("./",outputPath,sep="")
  }

  print("====================================================================")
  print("Ending ARTIVAnet procedure")
  print("====================================================================")

  print("====================================================================")
  print(paste("All results are saved in Folder:", outputPath ))
  print("====================================================================")
  
# End of the function ARTIVAnet()
  return(GlobalNetwork)
}
