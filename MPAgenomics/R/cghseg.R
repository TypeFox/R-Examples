### Model selection function with C constant
# @note used in bestSegmentationBM
modelSelection <- function(ERM, n, C)
{
  D <- length(ERM)
  pen <- C*(1:D)*(2.5+log(n/1:D))/n
  DHat <- which.min(ERM/n+pen)
  return(DHat)
}


# Method from Birge Massart to find an optimal segmentation of cghseg results.
#
# @param res output from cghseg segmentation (segmeanCO)
# @param n length of signal
#
# @return a list  containing
# \describe{
#   \item{D}{number of segments}
#   \item{m}{breakpoints}
#   \item{C}{penalty parameter}
#}
#
# @examples
# profil=c(rnorm(20,0,0.5),rnorm(10,-1,0.5),rnorm(40,1,0.5),rnorm(20,0,0.5))
# Kmax=10
# rescghseg <- cghseg:::segmeanCO(profil, Kmax=Kmax)
# bestSegmentationBM(rescghseg,length(profil))
#
bestSegmentationBM <- function(res,n)
{
  ERM=res$J.est
  Dmax = length(ERM)
  highD = (Dmax%/%2):Dmax
  model <- lm(ERM[highD]~highD)
  lambda <- -model$coefficient["highD"]
  DhatBest <- modelSelection(ERM, n ,lambda)
  mhatBest <- res$t.est[DhatBest, 1: DhatBest] 
  return(list(D = DhatBest, m =  mhatBest, C = lambda))
}

# Second Method from Birge Massart to find an optimal segmentation of cghseg results.
#
# @param lambdaList vector of constants in decreasing order
# @param resCGHSeg output from cghseg segmentation (segmeanCO)
# @param n length of signal
#
# @return a list  containing
# \describe{
#   \item{D}{number of segments}
#   \item{m}{breakpoints}
#   \item{C}{penalty parameter}
#}
#
# @examples
# profil=c(rnorm(20,0,0.3),rnorm(10,-1,0.3),rnorm(40,1,0.3),rnorm(20,0,0.3))
# Kmax=10
# rescghseg <- cghseg:::segmeanCO(profil, Kmax=Kmax)
# bestSegmentation(seq(10,0,length=1000),rescghseg,length(profil))
#
bestSegmentation <- function(lambdaList,resCGHSeg,n)
{
  Dhat <- sapply(lambdaList, modelSelection, ERM = resCGHSeg$J.est, n = n)
  maxJumpInDhat <- which.max(diff(Dhat))+1
  DhatBest <- modelSelection(resCGHSeg$J.est, n ,lambdaList[maxJumpInDhat]*2)
  mhatBest <- resCGHSeg$t.est[DhatBest, 1:DhatBest]
  return(list(D = DhatBest, m =  mhatBest, C = lambdaList[maxJumpInDhat]*2))
}


#
# This function launches the segmentation cghseg (from package cghseg).
# 
# @title segmentation function 
#
# @param signal a vector containing the signal.
# @param Kmax maximal number of segments
# @param position A vector containing the position of all elements of the signal (not necessary)
# @param plot if TRUE, plot the segmentation results
# @param verbose if TRUE print some informations
# 
# @return a list containing
# \describe{
#   \item{signal}{A vector containing the signal.}
#   \item{segmented}{A vector of the same size as signal containing the segmented values.}
#   \item{startPos}{The position of each probe.}
#   \item{segment}{A data.frame that summarizes the results of the segmentation. Each row is a different segment with the start position, end position, number of points in the signal and the value of the segment.}
# }
#
# @export
# 
# @author Quentin Grimonprez
# 
cghseg=function(signal,Kmax=10,position=NULL,plot=TRUE,verbose=TRUE)
{
  #signal
  if(missing(signal))
    stop("signal is missing.") 
  if(!is.numeric(signal) )#|| !is.vector(signal))
    stop("signal must be a vector of real.")
  
  #Kmax
  if(!is.numeric(Kmax) || (length(Kmax)>1))
  if(!is.wholenumber(Kmax))
    stop("Kmax must be a positive integer.")
  if(Kmax<=0)
    stop("Kmax must be a positive integer.")
  
  #position
  if(is.null(position))
    position=1:length(signal)
  if(!is.numeric(position) || !is.vector(position))
    stop("position must be a vector of real.")
  
  #order signal
  ord=order(position)
  position=position[ord]
  signal=signal[ord]  

  
  #doesn't tolerate NA values
  noNA=which(!is.na(signal))
  
  #Test if there is at least 2 valid points.
  if (length(noNA) < 2)
  {
    warning("Not enough point to segment signal")
    return(NULL)
  }
  
    
  #segmentation
  segmentation=cghseg:::segmeanCO(signal, Kmax=Kmax)
  
  #find best segmentation
  #best=bestSegmentationBM(segmentation,length(signal))
  lambdalist=seq(var(signal)*10,0,length=1000)
  best=bestSegmentation(lambdalist,segmentation,length(signal))
  
  cpt=c(0,best$m)
  means=unlist(lapply(1:length(best$m),FUN=function(i){return(mean(signal[(cpt[i]+1):cpt[i+1]]))}))
  
  if(plot)
  {    
    #plot data
    plot(position,signal,pch=".",xlab="Position",ylab="signal",ylim=c(0,6))
    
    #plot segments 
    for(i in 1:length(best$m))
      lines(c((position[noNA])[cpt[i]+1],(position[noNA])[cpt[i+1]]),rep(means[i],2),col="red",lwd=3)  
  }
  
  #create segmented signals
  nbPtsSeg=diff(cpt)
  segmentedSignal=rep(NA,length(signal))  
  segmentedSignal[noNA]=unlist(lapply(1:length(nbPtsSeg),FUN=function(i)
  {
    rep(means[i],nbPtsSeg[i])
  }))
  
  return(list(signal=as.matrix(signal),
              segmented=as.matrix(segmentedSignal),
              startPos=position,
              segment=data.frame(start=(position[noNA])[cpt[-length(cpt)]+1],
                                 end=(position[noNA])[cpt[-1]],
                                 points=nbPtsSeg,
                                 means=means)))
}



#########################################################################


#
# This function launches the segmentation cghseg.
#
# @title segmentation function 
#
# @param dataSetName The name of the data-set folder (it must correspond to a folder name in rawData folder.).
# @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param chromosome A vector with the chromosomes to be segmented. 
# @param Kmax maximum number of segments.
# @param onlySNP If TRUE, only the copy-number for SNPs positions will be returned (default=TRUE).
# @param listOfFiles A vector containing the names of the files in dataSetName folder for which the copy number profiles will be segmented (default is all the files).
# @param savePlot if TRUE, graphics of the segmented CN signal will be saved in the figures/dataSetName/segmentation/CN folder. (default=TRUE).
# @param verbose if TRUE print some informations
# 
# @return a list containing
# \describe{
#   \item{copynumber}{A vector containing  the copynumber signal.}
#   \item{segmented}{A vector of the same size as copynumber containing the segmented values.}
#   \item{startPos}{The position of each probes.}
#   \item{chromosome}{A vector of the same size as copynumber containing the chromosome number.}
#   \item{featureNames}{Names of the probes.}
#   \item{sampleNames}{The name of the signal.}
#   \item{segment}{A data.frame that summarizes the results of the segmentation. Each row is a different segment with the chromosome, start position, end position, number of probes in the signal and the value of the segment.}
# }
#
# @export
# 
# @author Quentin Grimonprez
# 
CGHSEGaroma=function(dataSetName,normalTumorArray,chromosome=1:22,Kmax=10,listOfFiles=NULL,onlySNP=TRUE,savePlot=TRUE,verbose=TRUE)
{
  require(R.devices)
  allpkg=TRUE
  if(!suppressPackageStartupMessages(require("aroma.affymetrix", quietly=TRUE) ) )
  {
    cat("Package not found: aroma.affymetrix. For download it:\n")
    cat("source(\"http://www.braju.com/R/hbLite.R\")\n")
    cat(" hbLite(\"sfit\")\n")
    cat("source(\"http://bioconductor.org/biocLite.R\")\n")
    cat("biocLite(\"affxparser\")\n")
    cat("biocLite(\"DNAcopy\")\n")
    cat("biocLite(\"aroma.light\")\n")
    #     cat("source(\"http://aroma-project.org/hbLite.R\")\n")
    cat("install.packages(\"aroma.affymetrix\")\n")
    allpkg=FALSE
  }

  
  if(!suppressPackageStartupMessages(require("aroma.cn", quietly=TRUE) ) )
  {
    cat("Package not found: aroma.cn. For download it:\n")
    cat("install.packages(\"aroma.cn\")\n") 
    allpkg=FALSE
  }
  
  if(!allpkg)
    stop("You have to install some packages : Follow the printed informations.")
  
  require(aroma.core)
  require(R.filesets)
  
  if(!("totalAndFracBData"%in%list.files()))
    stop("There is no \"totalAndFracBData\", check if you are in the good working directory or if you have run the signalPreProcess function before.")
  
  ###check the arguments
  #dataSetName
  if(missing(dataSetName))
    stop("dataSetName is missing.")
  if(!is.character(dataSetName))
    stop("dataSetName must be the name of a folder in \"rawData\" folder.")
  if(!(dataSetName%in%list.files("rawData")))
    stop("dataSetName must be the name of a folder in \"rawData\" folder.")
  #onlySNP
  if(!is.logical(onlySNP))
    stop("onlySNP must be a boolean.")
  
  #Kmax
  if(!is.numeric(Kmax) || (length(Kmax)>1))
    if(!is.wholenumber(Kmax))
      stop("Kmax must be a positive integer.")
  if(Kmax<=0)
    stop("Kmax must be a positive integer.")
  
  #check if we are in a normal-tumor study or in a single array study
  singleStudy=TRUE
  if(!missing(normalTumorArray))
    singleStudy=FALSE
  
  #################### import the dataSet to have the name of all the files
  
  #path where find the CN data
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  #load CN
  dsC <- aroma.core::AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
  
  ################### check normalTumorArray
  if(!singleStudy)
  {
    #normalTumorArray
    if(is.character(normalTumorArray))
      normalTumorArray=read.csv(normalTumorArray)
    else
    {
      if(!is.data.frame(normalTumorArray))
        stop("normalTumorArray must be either the path to the normalTumorArray csv file or a data.frame containing the data.\n")
    }
    
    #check normalTumorArray
    if(!("normal"%in%colnames(normalTumorArray)) || !("tumor"%in%colnames(normalTumorArray)))
      stop("normalTumorArray doesn't contain a column \"normal\" or \"tumor\".\n")
    
    #     isArrayComplete=sapply(getNames(dsC),FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
    #     if(sum(isArrayComplete)!=length(isArrayComplete))
    #       stop("normalTumorArray doesn't contain all the filenames of dataSetName.")
  }
  
  ###### check listOfFiles
  pos=c()
  if(is.null(listOfFiles) || missing(listOfFiles))
  {
    listOfFiles=R.filesets::getNames(dsC)
    pos=1:length(dsC)
  }
  else
  {
    #check the format
    if(!is.vector(listOfFiles) || !is.character(listOfFiles))
      stop("listOfFiles must be a vector of string.")
    listOfFiles=unique(listOfFiles)
    #check if all the files of listOfFiles are in the folder
    pos=sapply(listOfFiles,match,R.filesets::getNames(dsC))#position of the files of listOfFiles in the folder
    if(length(which(pos>0))!=length(pos))
      stop("Wrong name of files in listOfFiles")
  }  
  
  #check if file of listOfFiles are in normalTumorArray
  if(!singleStudy)
  {
    isArrayComplete=sapply(listOfFiles,FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
    if(sum(isArrayComplete)!=length(isArrayComplete))
      stop("normalTumorArray doesn't contain all the filenames you specified in listOfFiles parameter.")
  }
  
  #if single array study, we just reduce dsC to pos  
  if(!singleStudy)
  {      
    #if normal-tumor study, we need the tumor and normal files
    
    #we obtain the complementary files  
    compFiles=getComplementaryFile(listOfFiles,normalTumorArray)
    allFiles=unique(c(listOfFiles,compFiles))
    
    #index of the files
    pos=sapply(allFiles,FUN=function(x,dsC){which(R.filesets::getNames(dsC)==x)},dsC)
    tag=getStatus(allFiles,normalTumorArray)
    
    pos=pos[which(tag=="tumor")]
  }
    
  ######################### END CHECK PARAMETERS
  
  ###bed files output
  #results are stored in the segmentation folder
  #check the existence of the segmentation folder
  if(!("segmentation"%in%list.files()))
    dir.create("segmentation");
  #check the existence of th dataSet folder in segmentation folder
  if(!(dataSetName%in%list.files("segmentation")))
    dir.create(paste0("segmentation/",dataSetName));
  if(!("CN"%in%list.files(paste0("segmentation/",dataSetName))))
    dir.create(paste0("segmentation/",dataSetName,"/CN"));
  
  figPath <- Arguments$getWritablePath(paste0("figures/",dataSetName,"/segmentation/CN/"));
  
  #names of the files to segment
  names=R.filesets::getNames(dsC)[pos]
  
  output=lapply(names,FUN=function(name)
  {
    cghArg=list()
    for(chr in chromosome)
    {
      #get the copy-number for 1 chr
      if(!singleStudy)
        CN=getCopyNumberSignal(dataSetName,chr,normalTumorArray,onlySNP,name,verbose=FALSE)
      else
        CN=getCopyNumberSignal(dataSetName,chr,onlySNP=onlySNP,listOfFiles=name,verbose=FALSE)
      
      CN=CN[[paste0("chr",chr)]]
      
      gc()
      
      #segmentation
      if (length(which(!is.na(as.vector(CN[,3]))))<2)
      {
        if (chr==24)
        {
          cat(paste0("Cannot segment file ",name," chromosome Y (24) : gender = XX\n"))
        } else {
          cat(paste0("Cannot segment file ",name," chromosome ",chr,  ": less than 2 points in the signal\n"))
        }
      } else {
        
        cat(paste0("Segmentation of file ",name," chromosome ",chr,"..."))
        seg=cghseg(as.vector(CN[,3]),Kmax,CN$position,plot=savePlot,verbose=FALSE)
        cat("OK\n")
        
        if(savePlot)
        {
          figName <- sprintf("%s,%s", name, chr);
          pathname <- filePath(figPath, sprintf("%s.png", figName));
          width <- 1280;
          aspect <- 0.6*1/3;
          fig <- R.devices::devNew("png", pathname, label=figName, width=width, height=2*aspect*width);
          plot(NA,xlim=c(min(CN$position),max(CN$position)), ylim=c(0,6),xlab="Position", main=figName,ylab="CN", pch=".")
          points(CN$position, seg$signal, pch=".");
          for(i in 1:nrow(seg$segment))
            lines(c(seg$segment$start[i],seg$segment$end[i]),rep(seg$segment$means[i],2),col="red",lwd=3)
          R.devices::devDone();
        }
        
        
        
        #concatenate the results
        cghArg=list(copynumber=rbind(cghArg$copynumber,seg$signal),
                    segmented=rbind(cghArg$segmented,seg$segmented),
                    startPos=c(cghArg$startPos,CN$position),
                    chromosome=c(cghArg$chromosome,CN$chromosome),
                    sampleNames=name,
                    featureNames=c(cghArg$featureNames,CN$featureNames),
                    segment=data.frame(chrom=c(as.character(cghArg$segment$chrom),rep(paste0("chr",CN$chromosome[1]),length(seg$segment$start))),
                                       chromStart=c(cghArg$segment$chromStart,seg$segment$start),
                                       chromEnd=c(cghArg$segment$chromEnd,seg$segment$end),
                                       probes=c(cghArg$segment$probes,seg$segment$points),
                                       means=c(cghArg$segment$means,seg$segment$means)))
      }
    }
    
    #write in .bed
    write.table(cghArg$segment,file=paste0("segmentation/",dataSetName,"/CN/",name,",segmentation.bed"),sep="\t",row.names=FALSE)
    
    return(cghArg)
  })
  
  names(output)=names  
  return(output)            
}

######################################
#'
#' This function launches the segmentation of a signal.
#'
#' @title segmentation function 
#'
#' @param signal a vector containing the signal.
#' @param method method of segmentation, either "PELT" or "cghseg".
#' @param Rho For method="PELT", vector containing all the penalization values to test for the segmentation. If no values are provided, default values will be used.
#' @param Kmax For method="cghseg", maximal number of segments.
#' @param position A vector containing the position of all elements of the signal (not necessary)
#' @param plot if TRUE, plot the segmentation results
#' @param verbose if TRUE print some informations
#' 
#' @return a list containing
#' \describe{
#'   \item{signal}{A vector containing the signal.}
#'   \item{segmented}{A vector of the same size as signal containing the segmented values.}
#'   \item{startPos}{The position of each probe.}
#'   \item{segment}{A data.frame that summarizes the results of the segmentation. Each row is a different segment with the start position, end position, number of points in the signal and the value of the segment.}
#' }
#'
#' @export
#' 
#' @author Quentin Grimonprez
#' 
segmentation=function(signal,method=c("cghseg","PELT"),Rho=NULL,Kmax=10,position=NULL,plot=TRUE,verbose=TRUE)
{
  method <- match.arg(method)
  seg=switch(method,
             PELT=PELT(signal,Rho,position,plot,verbose),
             cghseg=cghseg(signal,Kmax,position,plot,verbose))
  
  return(seg)
}

######################################
#'
#' This function launches the segmentation process using the aroma architecture.
#'
#' @title segmentation function 
#'
#' @param dataSetName The name of the data-set folder (it must correspond to a folder name in rawData folder.).
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param chromosome A vector with the chromosomes to be segmented.
#' @param method method of segmentation, either "PELT" or "cghseg".
#' @param Rho For method="PELT", vector containing all the penalization values to test for the segmentation. If no values are provided, default values will be used.
#' @param Kmax For method="cghseg", maximal number of segments.
#' @param onlySNP If TRUE, only the copy-number for SNPs positions will be returned (default=TRUE).
#' @param listOfFiles A vector containing the names of the files in dataSetName folder for which the copy number profiles will be segmented (default is all the files).
#' @param savePlot if TRUE, graphics of the segmented CN signal will be saved in the figures/dataSetName/segmentation/CN folder. (default=TRUE).
#' @param verbose if TRUE print some informations
#' 
#' @return a list containing
#' \describe{
#'   \item{copynumber}{A vector containing  the copynumber signal.}
#'   \item{segmented}{A vector of the same size as copynumber containing the segmented values.}
#'   \item{startPos}{The position of each probes.}
#'   \item{chromosome}{A vector of the same size as copynumber containing the chromosome number.}
#'   \item{featureNames}{Names of the probes.}
#'   \item{sampleNames}{The name of the signal.}
#'   \item{segment}{A data.frame that summarizes the results of the segmentation. Each row is a different segment with the chromosome, start position, end position, number of probes in the signal and the value of the segment.}
#' }
#'
#' @export
#' 
#' @author Quentin Grimonprez
#' 
segmentationAroma=function(dataSetName,normalTumorArray,chromosome=1:22,method=c("cghseg","PELT"),Kmax,Rho=NULL,listOfFiles=NULL,onlySNP=TRUE,savePlot=TRUE,verbose=TRUE)
{
  method <- match.arg(method)
  seg=switch(method,
             PELT=PELTaroma(dataSetName,normalTumorArray,chromosome,Rho,listOfFiles,onlySNP,savePlot,verbose),
             cghseg=CGHSEGaroma(dataSetName,normalTumorArray,chromosome,Kmax,listOfFiles,onlySNP,savePlot,verbose))
  
  return(seg)
}