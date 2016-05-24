#'
#' This function applies the PELT method to segment each signal of the dataset and launches CGHcall for calling segments and detect aberrations.
#' Results will be stored in a text file in the segmentation folder of the aroma architecture.
#'
#' @title Segment a copy-number signal and call the found segments.
#'
#' @param dataSetName name of the data-set folder in the rawData folder containing the signals to use.
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files.
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param chromosome A vector containing the chromosomes to segment.
#' @param method method of segmentation, either "PELT" or "cghseg".
#' @param Rho For method="PELT", vector containing all the penalization values to test for the segmentation. If no values are provided, default values will be used.
#' @param Kmax For method="cghseg", maximal number of segments.
#' @param listOfFiles A vector containing the names of the files from the dataSetName to use.
#' @param onlySNP If TRUE, only the SNP probes will be used.
#' @param savePlot If TRUE, save the segmented signal in figures folder.
#' @param nclass The number of levels to be used for calling. Either 3 (loss, normal, gain), 4 (including amplifications), 5 (including double deletions) (default=3).
#' @param cellularity Percentage of tumored cells in the sample (default=1).
#' @param ... Other parameters of CGHcall function
#'
#' @return a data.frame containg columns :
#' \describe{
#'   \item{sampleNames}{Name of the file.}
#'   \item{chrom}{The chromosome of the segment.}
#'   \item{chromStart}{The starting position (in bp) of a segment. This position is not included in the segment.}
#'   \item{chromEnd}{The ending position (in bp) of a segment.This position is included in the segment.}
#'   \item{probes}{Number of probes in the segment.}
#'   \item{means}{Mean of the segment.}
#'   \item{calls}{The calling of segment ("double loss", "loss", "normal", "gain" or "amplification").}
#' }
#'
#' @examples
#' #DO NOT EXECUTE before reading the vignette
#' # seg1=cnSegCallingProcess("data1",normalTumorArray,chromosome=20:21)
#' # seg2=cnSegCallingProcess("data2",chromosome=20:21)
#'
#' @export
#' 
#' @author Quentin Grimonprez
#' 
cnSegCallingProcess=function(dataSetName,normalTumorArray,chromosome=1:22,method=c("cghseg","PELT"),Rho=NULL,Kmax=10,listOfFiles=NULL,onlySNP=TRUE,savePlot=TRUE,nclass=3,cellularity=1,...)
{
  require(R.devices)
  method <- match.arg(method)
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
  #   else
  #     cat("Package aroma.affymetrix loaded.\n")
  
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
    stop("onlyDNP must be a boolean.")

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
    
#     isArrayComplete=sapply(R.filesets::getNames(dsC),FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
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
  
  #Rho
  if(missing(Rho) || is.null(Rho))
    Rho=c(seq(0.1,2,by=0.1),seq(2.2,5,by=0.2),seq(5.5,10,by=0.2),seq(11,16,by=1),seq(18,36,by=2),seq(40,80,by=4))

  ######################### END CHECK PARAMETERS
  
  figPath <- Arguments$getWritablePath(paste0("figures/",dataSetName,"/segmentation/CN/"));
  
  #names of the files to segment
  names=R.filesets::getNames(dsC)[pos]
  
  output=lapply(names,FUN=function(name)
  {
    cghArg=list()
    segmentList=data.frame()
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

        seg=switch(method,
          PELT=PELT(as.vector(CN[,3]),Rho,CN$position,plot=TRUE,verbose=FALSE),
          cghseg=cghseg(as.vector(CN[,3]),Kmax,CN$position,plot=TRUE,verbose=FALSE))
          
        cat("OK\n")
        
        if(savePlot)
        {
          figName <- sprintf("%s,%s", name, chr);
          pathname <- filePath(figPath, sprintf("%s.png", figName));
          width <- 1280;
          aspect <- 0.6*1/3;
          fig <- R.devices::devNew("png", pathname, label=figName, width=width, height=2*aspect*width);
          plot(NA,xlim=c(min(CN$position),max(CN$position)), ylim=c(0,6),xlab="Position", main=figName,ylab="CN", pch=".")
          points(CN$position, CN[,3], pch=".");
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
                    featureNames=c(cghArg$featureNames,CN$featureNames))
  
        #output for the bed file      
        segmentList=data.frame(chrom=c(as.character(segmentList$chrom),rep(paste0("chr",CN$chromosome[1]),length(seg$segment$start))),
                          chromStart=c(segmentList$chromStart,seg$segment$start),
                          chromEnd=c(segmentList$chromEnd,seg$segment$end),
                          probes=c(segmentList$probes,seg$segment$points),
                          means=c(segmentList$means,seg$segment$means))
      }
    }
        
    if (length(cghArg)>0)
    {
    #launch calling process
    cghArg=callingProcess(cghArg,nclass=nclass,cellularity=cellularity,verbose=FALSE,...)
    
    #get the calling for each segment
    calls=apply(segmentList,1,FUN=function(segment)
      {

        #we restrict to the chromosome of the segment
        indChr=which(cghArg$chromosome==as.numeric(substr(as.character(segment[1]),4,nchar(as.character(segment[1])))))

        #we search the starting position of the segment
        indPos=which(cghArg$startPos[indChr]==(as.numeric(segment[2])))[1]
      
        #print(indPos)
        # calling of the start of the segment=calling of the segment
        return(cghArg$calls[indChr[indPos]]) 
      })
    
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
    #output=segment+calls
    output=data.frame(segmentList,calls=calls)
    #write in .bed
    write.table(output,file=paste0("segmentation/",dataSetName,"/CN/",name,",calls.bed"),sep="\t",row.names=FALSE)
    
    return(data.frame(sampleNames=rep(name,nrow(output)),output))
   }
      
  })

  return(do.call(rbind,output))
}

##########################################################################################################################################

#'
#' This function filters the output of a segmentation and label process.
#' It allows to keep only segments over a minimal length or containing at least a minimal number of probes.
#'
#' @title Filter segments
#'
#' @param segmentList A data.frame containing a description of segments, it must have at least columns named "chromStart", "chromEnd", "probes" and "calls".
#' (see the output of \link{cnSegCallingProcess} function).
#' @param minLength The minimum length (in bp) for a segment. All the shorter segments are removed.
#' @param minProbes The minimum number of probes for a segment. All the segments with less probes are removed.
#' @param keptLabel Vector of labels to keep. Only segment with one of the specified label will be kept.
#'
#' @return a data.frame of the same format as segmentList.
#'
#' @export
#' 
#' @author Quentin Grimonprez
#' 
filterSeg=function(segmentList,minLength=1,minProbes=1,keptLabel=c("loss","gain"))
{
  ########check parameters
  #segmentList
  if(missing(segmentList))
    stop("segmentList is missing.")
  if(!is.data.frame(segmentList))
    stop("segmentList must be a data.frame.")
  if(!("calls"%in%names(segmentList)))
    stop("The column \"calls\" from segmentList is missing.")  
  if(!("probes"%in%names(segmentList)))
    stop("The column \"probes\" from segmentList is missing.")  
  if(!("chromStart"%in%names(segmentList)))
    stop("The column \"chromStart\" from segmentList is missing.")  
  if(!("chromEnd"%in%names(segmentList)))
    stop("The column \"chromEnd\" from segmentList is missing.")  
  #minLength
  if(!is.numeric(minLength))
    stop("minLength must be a positive integer.")
  if(minLength<0)
    stop("minLength must be a positive integer.")
  if(!is.wholenumber(minLength))
    stop("minLength must be a positive integer.")
  #minProbes
  if(!is.numeric(minProbes))
    stop("minProbes must be a positive integer.")
  if(minProbes<0)
    stop("minProbes must be a positive integer.")
  if(!is.wholenumber(minProbes))
    stop("minProbes must be a positive integer.")
  #keptLabel
  if(!is.character(keptLabel) || !is.vector(keptLabel))
    stop("keptLabel must be a vector of chracter.")
  
  #########
  
  #keep only segments with the specified label
  ind=c()
  for(label in keptLabel)
    ind=c(ind,which(segmentList$calls==label))
  segmentList=segmentList[ind,]  
  #remove the segments shorter than minLength pb
  segmentList=segmentList[segmentList$chromEnd-segmentList$chromStart>=minLength,]
  #remove the segments containing not enough probes
  segmentList=segmentList[segmentList$probes>=minProbes,]
  
  return(segmentList)
}