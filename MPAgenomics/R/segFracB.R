#'
#' This function launches the segmentation of allele B fraction only for heterozygous SNPs.
#'
#' @title segmentation function for the allele B fraction
#'
#' @param dataSetName The name of the data-set folder (it must correspond to a folder name in rawData folder.).
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files.
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param chromosome  A vector with the chromosomes to be segmented. 
#' @param method method of segmentation, either "PELT" or "cghseg".
#' @param Rho For method="PELT", vector containing all the penalization values to test for the segmentation. If no values are provided, default values will be used.
#' @param Kmax For method="cghseg", maximal number of segments.
#' @param listOfFiles A vector containing the names of the files in dataSetName folder for which the allele B profile is segmented (default is all the files).
#' @param savePlot if TRUE, graphics of the segmented allele B profile will be saved in the figures/dataSetName/segmentation/fracB folder. (default=TRUE).
#' @param verbose if TRUE print some informations
#' 
#' @return a data.frame where each row correspond to a different segment with columns :
#' \describe{
#'   \item{sampleNames}{The name of the signal.}
#'   \item{chromosome}{A vector of the same size as copynumber containing the chromosome number.}
#'   \item{chromStart}{The starting position of a segment.}
#'   \item{chromEnd}{The ending position of a segment.}
#'   \item{probes}{The number of probes in the segment.}
#'   \item{means}{Means of the segment.}
#' }
#'
#' @export
#' 
#' @author Quentin Grimonprez
#' 
segFracBSignal=function(dataSetName,normalTumorArray,chromosome=1:22,method=c("cghseg","PELT"),Rho=NULL,Kmax=10,listOfFiles=NULL,savePlot=TRUE,verbose=TRUE)
{
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

  if(!suppressPackageStartupMessages(require("aroma.cn", quietly=TRUE) ) )
  {
    cat("Package not found: aroma.cn. For download it:\n")
    cat("install.packages(\"aroma.cn\")\n") 
    allpkg=FALSE
  }
  
  if(!allpkg)
    stop("You have to install some packages : Follow the printed informations.")
  
  require(aroma.core)
  require(R.devices)
  require(R.filesets)
  
  if(!("totalAndFracBData"%in%list.files()))
    stop("There is no \"totalAndFracBData\", check if you are in the good working directory or if you have run the signalPreProcess function before.")
  if(!("callData"%in%list.files()))
    stop("There is no \"totalAndFracBData\", check if you are in the good working directory or if you have run the signalPreProcess function before.")
  
  ###check the arguments
  #dataSetName
  if(missing(dataSetName))
    stop("dataSetName is missing.")
  if(!is.character(dataSetName))
    stop("dataSetName must be the name of a folder in \"rawData\" folder.")
  if(!(dataSetName%in%list.files("rawData")))
    stop("dataSetName must be the name of a folder in \"rawData\" folder.")

  #################### import the dataSet to have the name of all the files
  #check if we are in a normal-tumor study or in a single array study
  #singleStudy=TRUE
  if(missing(normalTumorArray))
    stop("No normalTumorArray specified.\n Youd need to specify a normalTumorArray to extract allele B fraction")
  
  #path where find the CN data
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  #load CN
  dsC <- aroma.core::AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
    
  ################### check normalTumorArray
#   if(!singleStudy)
#   {
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
#       warning("normalTumorArray doesn't contain all the filenames of dataSetName.")
#   }
  
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
    
    if(!missing(normalTumorArray))
    {
      isArrayComplete=sapply(listOfFiles,FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
      if(sum(isArrayComplete)!=length(isArrayComplete))
        stop("normalTumorArray doesn't contain all the filenames you specified in listOfFiles parameter.")
    }

  }  
  
  #if single array study, we just reduce dsC to pos  
#   if(!singleStudy)
#   {  
#     #if normal-tumor study, we need the tumor and normal files
#     
    #we obtain the complementary files  
    compFiles=getComplementaryFile(listOfFiles,normalTumorArray)
    allFiles=unique(c(listOfFiles,compFiles))
    
    #index of the files
    pos=sapply(allFiles,FUN=function(x,dsC){which(R.filesets::getNames(dsC)==x)},dsC)
    tag=getStatus(allFiles,normalTumorArray)
    
    pos=pos[which(tag=="normal")]

    if (length(pos)!=length(which(tag=="tumor")))
    {
      stop("normalTumorArray must contain one unique normal sample per tumor sample to extract allele B fraction.")
    }

    
#   }
 
  #Rho
  if(missing(Rho) || is.null(Rho))
    Rho=c(seq(0.1,2,by=0.1),seq(2.2,5,by=0.2),seq(5.5,10,by=0.2),seq(11,16,by=1),seq(18,36,by=2),seq(40,80,by=4))
  
  ######################### END CHECK PARAMETERS
  ###bed files output
  #results are stored in the segmentation folder
  #check the existence of the segmentation folder
  if(!("segmentation"%in%list.files()))
    dir.create("segmentation");
  #check the existence of th dataSet folder in segmentation folder
  if(!(dataSetName%in%list.files("segmentation")))
    dir.create(paste0("segmentation/",dataSetName));
  
  figPath <- Arguments$getWritablePath(paste0("figures/",dataSetName,"/segmentation/fracB"));
  
  #names of the files to segment
  names=R.filesets::getNames(dsC)[pos]
  
  output=lapply(names,FUN=function(name)
  {
    segment=data.frame()
    for(chr in chromosome)
    {
      #get the fracB for 1 chr
     # if(!singleStudy)
        fracB=getFracBSignal(dataSetName,chromosome=chr,normalTumorArray,listOfFiles=name,verbose=verbose)
      #else
      #  fracB=getFracBSignal(dataSetName,chromosome=chr,listOfFiles=name,verbose=verbose)

      fracB=fracB[[paste0("chr",chr)]]$tumor
      gc()
      #get the genotype calls for 1 chr
      geno=getGenotypeCalls(dataSetName,chromosome=chr,listOfFiles=name,verbose=verbose)
      geno=geno[[paste0("chr",chr)]]
      gc()
      
      ind=which(geno[,3]=="AB")

      rm(geno)
      gc()
      
      fracB=fracB[ind,]

      fracB[,3]=symmetrizeFracB(fracB[,3])

      
      #segmentation
      cat(paste0("Segmentation of file ",name," chromosome ",chr,"..."))
      
      if (is.null(fracB[,3]) || length(fracB[,3])<2){
        cat("to few point to segment \n")
      } else {
        seg=switch(method,
          PELT=PELT(fracB[,3],Rho,position=fracB$position,plot=TRUE,verbose=verbose),
          cghseg=cghseg(fracB[,3],Kmax,position=fracB$position,plot=TRUE,verbose=verbose))
        cat("OK\n")
        
        if(savePlot)
        {
          figName <- sprintf("%s,%s,%s", name, "fracB,chr",chr);
          pathname <- filePath(figPath, sprintf("%s.png", figName));
          width <- 1280;
          aspect <- 0.6*1/3;
          fig <- R.devices::devNew("png", pathname, label=figName, width=width, height=2*aspect*width);
          plot(NA,xlim=c(min(fracB$position),max(fracB$position)), ylim=c(0,1),xlab="Position", main=figName,ylab="Allele B fraction", pch=".")
          points(fracB$position, fracB[,3], pch=".");
          for(i in 1:nrow(seg$segment))
            lines(c(seg$segment$start[i],seg$segment$end[i]),rep(seg$segment$means[i],2),col="red",lwd=3)
          R.devices::devDone();
        }
      
      
      
      #concatenate the results
#       cghArg=list(fracB=rbind(cghArg$copynumber,seg$signal),
#                   segmented=rbind(cghArg$segmented,seg$segmented),
#                   startPos=c(cghArg$startPos,fracB$position),
#                   chromosome=c(cghArg$chromosome,rep(fracB$chromosome,length(seg$signal))),
#                   sampleNames=name,
#                   featureNames=c(cghArg$featureNames,fracB$featureNames),
#                   segment=data.frame(chrom=c(as.character(cghArg$segment$chrom),unlist(lapply(rep(fracB$chromosome,length(seg$segment$means)),as.character))),
#                                      chromStart=c(cghArg$segment$chromStart,seg$segment$start),
#                                      chromEnd=c(cghArg$segment$chromEnd,seg$segment$end),
#                                      probes=c(cghArg$segment$probes,seg$segment$points),
#                                      means=c(cghArg$segment$means,seg$segment$means)))
      
        segment=data.frame(chrom=c(as.character(segment$chrom),rep(paste0("chr",chr),length(seg$segment$start))),
                           chromStart=c(segment$chromStart,seg$segment$start),
                           chromEnd=c(segment$chromEnd,seg$segment$end),
                           probes=c(segment$probes,seg$segment$points),
                           means=c(segment$means,seg$segment$means))
      }
    }
    
    
    #write in .bed
    write.table(segment,file=paste0("segmentation/",dataSetName,"/",name,",symFracB,segmentation.bed"),sep="\t",row.names=FALSE)
    segment=cbind(rep(name,nrow(segment)),segment)
    names(segment)[1]="sampleNames"
    return(segment)
  })
  
  output=do.call(rbind,output)
  row.names(output)=NULL
  
  return(output)            
}
