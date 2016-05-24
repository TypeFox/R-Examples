ARTIVAsubnetAnalysis <-
function( ARTIVAsubnet=NULL, CPpostDist=NULL,CPsamples=NULL,coefSamples=NULL,TFnumber=NULL, segMinLength=2, edgesThreshold=0.5, burn_in=NULL, CPpos=NULL,targetData=NULL, parentData=NULL,targetName=NULL,parentNames=NULL, savePictures=TRUE,saveEstimations=TRUE,outputPath=NULL,layout="fruchterman.reingold", silent=FALSE, inARTIVAsubnet=FALSE , onepage= FALSE){
   
  if(!(is.null(ARTIVAsubnet))){
    TFnumber=ARTIVAsubnet$GLOBvar$q
    targetData=ARTIVAsubnet$targetData
    parentData=ARTIVAsubnet$parentData
    if(is.null(targetName))targetName=ARTIVAsubnet$GLOBvar$targetName
    if(is.null(parentNames))parentNames=ARTIVAsubnet$GLOBvar$parentNames
    if(is.null(burn_in))burn_in=ARTIVAsubnet$OUTvar$burn_in
    niter=ARTIVAsubnet$GLOBvar$niter
    CPpostDist=ARTIVAsubnet$CPpostDist
    CPsamples=ARTIVAsubnet$Samples$CP
    coefSamples=ARTIVAsubnet$Samples$coeff
  }else{
    if(is.null(CPpostDist) ){
      stop("Please specify either parameter CPpostDist or parameter ARTIVAsubnet ")
    }
    if(is.null(CPsamples) ){
      stop("Please specify either parameter CPsamples or parameter ARTIVAsubnet ")
    }
    if(is.null(coefSamples) ){
       stop("Please specify either parameter coefSamples or parameter ARTIVAsubnet ")
   }
    if(is.null(TFnumber) ){
       stop("Please specify either parameter TFnumber or parameter ARTIVAsubnet ")
   }
    niter=dim(CPsamples)[1]
    if(is.null(burn_in))burn_in=round(niter/4)

  }
  
  CPsamples=CPsamples[burn_in:niter,]
  coefSamples=coefSamples[burn_in:niter,]


   ## CP position
  if(is.null(CPpos)){
    CPpos=CPposition(nbCP=ARTIVAsubnet$CPpostDist$estimatedCPnumber,CPpositionPostDist=CPpostDist$CPposition,CPsamples=CPsamples,segMinLength=segMinLength)
  }else{
    if( sum(((CPpos[-1]-CPpos[-length(CPpos)])>=segMinLength)==FALSE)>0){
      stop("Please specify compatible CP position vector (CPpos) and minimum temporal segment length(segMinLength)")
    }
  }

    ## General information concerning the analyzed data set
  if(!inARTIVAsubnet & !silent)
  {
    print("========================")
    print("ARTIVA subnet Analysis")
    print(paste(" *** Name of the analyzed target gene:", targetName))
    print(paste(" **** Minimal length to define a temporal segment:", segMinLength, "time points"))
    print("========================")
  }

  
  
  if(savePictures | saveEstimations ){
    if(is.null(outputPath)){
      # Output directory path (if it doesn t exist create it)
      path =getwd()
      if(.Platform$OS.type == "unix"){
        if(! "ARTIVAsubnet" %in% system("ls" ,intern=TRUE)) {
        system(paste("mkdir ", path,"/ARTIVAsubnet", sep="")) }
      }else{# if(.Platform$OS.type == "unix"){
        shell("mkdir ARTIVAsubnet", intern = TRUE,mustWork =NA)
      }
      outputPath="ARTIVAsubnet/"
    }else{
      if(.Platform$OS.type == "unix"){
        if(! outputPath  %in% system( "ls ",intern=TRUE)) {
          system(paste("mkdir ",outputPath, sep=""))
        }
      }else{# if(.Platform$OS.type == "unix"){
        shell(paste("mkdir ",outputPath, sep=""), intern = TRUE,mustWork =NA)
      }
    }
   }

############################## 
  ##  number of temporal segments 
  nbSegs=length(CPpos)-1
  #nbSegs=order(CPpostDist$CPnumber,decreasing=TRUE)[1]
  SegmentPostDist=NULL  
  TraceNetwork=NULL
  
  if(!silent){
    if(nbSegs==1){
      print("Only 1 temporal segment was identified (maximizing the CPs posterior distribution)")
   }else{
      print(paste(nbSegs, "different temporal segments were identified (maximizing the CPs posterior distribution)"))
    }
  }
  if (!silent){
    if(inARTIVAsubnet){
      print("STEP 4: Computing the posterior distribution for the regulatory models in each temporal phase")
    }else{
      print("Computing the posterior distribution for the regulatory models in each temporal phase")
    }
  }
  
  ## edges analysis: compute the Edges posterior distribution within the selected segments, that is CP number with highest posterior probability and corresponding number of CP positions with highest posterior probability.
  SegmentPostDist=segmentModel.postDist(CPnumberPostDist=CPpostDist$CPnumber,CPpositionPostDist=CPpostDist$CPposition,TFnumber=TFnumber,TFnames=parentNames,CPsamples=CPsamples,coefSamples=coefSamples,segMinLength=segMinLength, edgesThreshold=edgesThreshold, CPpos=CPpos)
 
  if(sum(is.na(SegmentPostDist$edgesPostDist))>0){
    print("--- !!! WARNING MESSAGE !!! ----")
    print("Segment analysis was not possible for SOME segments because a lack of convergence, please use ARTIVAnet or ARTIVAsubnet again with more iterations")
    print("---------------------------------")
  }
 
  if(saveEstimations){
    if(.Platform$OS.type == "unix"){
      if(! "Estimations" %in% system(paste( "ls ",outputPath),intern=TRUE)) { system(paste("mkdir ",outputPath,"/Estimations",sep=""))}
    }else{# if(.Platform$OS.type == "unix"){
      tmp=getwd()
      setwd(outputPath)
      shell("mkdir Estimations", intern = TRUE,mustWork =NA)
      setwd(tmp)
    }
    outputResPath=paste(outputPath,"/Estimations/",sep="")
 
    write.table(as.matrix(SegmentPostDist$edgesPostDist),paste(outputResPath,"edgesPostDist_",targetName,".txt",sep=""))

    write.table(as.matrix(SegmentPostDist$edgesCoeff),paste(outputResPath,"edgesCoeff_",targetName,".txt",sep=""))
  }

    if (!silent){
      print(paste("Graphical representation of the ARTIVA results, analyzing the target gene", targetName, "and searching for potential regulatory interaction with", TFnumber, "parent genes"))
    }
 
  if(savePictures){
    if(.Platform$OS.type == "unix"){
      if(! "Pictures" %in% system(paste( "ls ",outputPath),intern=TRUE)) { system(paste("mkdir ",outputPath,"/Pictures",sep=""))}
    }else{# if(.Platform$OS.type == "unix"){
      tmp=getwd()
      setwd(outputPath)
      shell("mkdir Pictures", intern = TRUE,mustWork =NA)
      setwd(tmp)

    }
    
    outputPicturesPath=paste(outputPath,"/Pictures/",sep="")
  
    if(!silent)print(paste("Saving Pictures at", getwd(), "/", outputPicturesPath, sep=""))
    pdf(paste(outputPicturesPath,"ARTIVA_graphics_", targetName, ".pdf",sep=""), width = 10, height = 7)

  }else{
    nbpictures=3+2*nbSegs
    par(mfrow = c(2,ceiling(nbpictures/2)))
 }
 
  # Change the police size
  par(cex.main = 0.9,cex.lab = 0.8)

  if(!is.null(targetData) & !is.null(parentData)){
    traceGeneProfiles(targetData, parentData, onepage=savePictures)#, CPpos = c(CPstart, CPend))
  }
     
  if(savePictures)par(mfrow = c(1,2))
   
  plotCP.postDist(CPpostDist, targetName=targetName, onepage=FALSE,  estimatedCPpos=CPpos)
 
  if(savePictures){
    if(nbSegs == 1| nbSegs == 2){ par(mfrow = c(1,nbSegs))}
    if(nbSegs == 3){ par(mfrow = c(2,2))}
    if(nbSegs > 4){ par(mfrow = c(3,2))}
  }
  plotSegmentModel.postDist(SegmentPostDist$edgesPostDist,SegmentPostDist$CPpos, targetName=targetName, parentNames=parentNames, edgesThreshold=edgesThreshold, onepage=FALSE)
 
  ## Table with all information to trace the time varying regulatory networkss
  Target   = NULL
  CPstart    = NULL
  CPend  = NULL
  Parent   = NULL
  PostProb = NULL

  ## Information of the coefficients value for the interaction
  CoeffMean = NULL
  
  for(S in 1:nbSegs)
  {
    for(P in 1:length(parentNames))
    {
      Target = c(Target, targetName)
      CPstart = c(CPstart, SegmentPostDist$CPpos[S])
      CPend = c(CPend, SegmentPostDist$CPpos[S + 1] - 1)
      Parent = c(Parent, parentNames[P])
      PostProb = c(PostProb, SegmentPostDist$edgesPostDist[S, P])
      # Information of the coefficients values for the interaction
      CoeffMean = c(CoeffMean, SegmentPostDist$edgesCoeff[S, P])
    }
  }
 
  TraceNetwork = data.frame(Parent, Target, CPstart, CPend,  PostProb, CoeffMean, edgesThreshold, row.names = NULL)
 
  #par(mfrow = NULL)
  if(savePictures){par(mfrow = c(1,1))}
  #x11()
  traceNetworks(TraceNetwork, edgesThreshold, layout=layout, onepage=onepage )  
  geneNetworkSummary(TraceNetwork, edgesThreshold)

  if(savePictures){
    dev.off()
    #dev.off()
  }
  
  if(saveEstimations)
  {  
    write.table(TraceNetwork, file = paste(outputPath,"/Estimations/ARTIVA_SubNetwork_", targetName, ".txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  
  ## return all samples	
  return(list(nbSegs=nbSegs, CPposition=SegmentPostDist$CPpos, SegmentPostDist=SegmentPostDist, network = TraceNetwork))
}
