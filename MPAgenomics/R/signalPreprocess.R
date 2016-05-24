#' low-level normalization process for estimating raw copy-numbers and allele B fraction.
#' 
#' 
#' @title Normalization process
#' @param dataFolder Name of the data set.
#' @param chipType Type of the chip used for the data.
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files.
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param genotypeCallsMethod method used for genotypage, default is "naive".
#' @param savePlot If TRUE, graphics of the CN signal and allele B fraction signal will be saved in the figures folder.
#' @param tags Common tag which appears in the different file names (cdf, ugp, ufl) of the chip. For no tag, use tags=NULL (default = NULL). See details for more information.
#' @return NULL
#'
#' @details The aroma architecture must be respected:
#'    <working directory>
#'     +- annotationData/
#'     |  +- chipTypes/
#'     |     +- <chipType>/ <-- must match exactly the name of the CDF file (fullname minus tags)
#'     |        +- CDF file(s) and other annotation (possibly subdirectories)
#'     |
#'     +- rawData/
#'        +- <nameOfDataSet>/
#'           +- <chipType>/ <-- must match exactly a chip type folder under annotationData/
#'              +- CEL files
#' 
#' 
#' All the cdf chip file names must follow the following rule : <chipType>,<Tags>.cdf
#' 
#' Multiples tags must be separated by  a comma. If there is no tag, the pattern is <chipType>.cdf
#' 
#' @author Quentin Grimonprez
#' 
#' @export
SignalNormalization<-function(dataFolder,chipType,normalTumorArray,genotypeCallsMethod="naive",savePlot=TRUE, tags=NULL)
{
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
#   else
#     cat("Package aroma.cn loaded.\n")

  if(!allpkg)
    stop("You have to install some packages : Follow the printed informations.")
  
  ###check required arguments
  #check architecture
  if(!("rawData"%in%list.files()) || !("annotationData"%in%list.files()) )
    stop("You have to respect the aroma architecture. See the documentation.")
  
  #dataFolder
  if(missing(dataFolder))
    stop("dataFolder is missing.")
  if(!is.character(dataFolder))
    stop("dataFolder must be of type character.")
  if(!(dataFolder%in%list.files("rawData")))
    stop(paste0(dataFolder," not found in rawData Folder."))
  
  #chipType
  if(missing(chipType))
    stop("chipType is missing.")
  if(!is.character(chipType))
    stop("chipType must be of type character.")
  if(!(chipType%in%list.files("annotationData/chipTypes")))
    stop(paste0(chipType," not found in annotationData/chipTypes Folder.")) 

  #check if ugp or ufl files with the specified tags exist
  file=list.files(paste0("annotationData/chipTypes/",chipType))
  if(!is.null(tags))
  {
    if(!is.character(tags))
      stop("tags must be a string.")
    
    tagfile=file[grep(tags,file)]
    
    tagfile=tolower(tagfile)
    cdffile=tagfile[grep(".cdf",tagfile)]
    ugpfile=tagfile[grep(".ugp",tagfile)]
    uflfile=tagfile[grep(".ufl",tagfile)]
    if(length(cdffile)==0)
      stop(paste0("No cdf files with tag ",tags))
    if(length(ugpfile)==0)
      stop(paste0("No ugp files with tag ",tags))
    if(length(uflfile)==0)
      stop(paste0("No ufl files with tag ",tags))
  }
  
  #genotypeCallsMethod
  if(genotypeCallsMethod!="naive")
    stop("genotypeCallsMethod must be \"naive\".")
  
  #check if we are in a normal-tumor study or in a single array study
  singleStudy=TRUE
  if(missing(normalTumorArray))
    cat("No normalTumorArray specified.\n A normalization and a genotype calls will be processed on all the data.\n")
  else
  {
    cat("A normal-tumor normalization will be processed with genotypage calls only for normal case and a tumorboost normalization for the tumor allele B fraction signal.")
    singleStudy=FALSE
  }

  ##check normalTumorArray
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
    
    #check is the file contains all the file
    fileWithoutCELExtension=sapply(list.files(paste0("rawData/",dataFolder,"/",chipType)),FUN=function(x){x=gsub(".cel","",x);x=gsub(".CEL","",x)})
    isArrayComplete=sapply(fileWithoutCELExtension,FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
    if(sum(isArrayComplete)!=length(isArrayComplete))
      stop("normalTumorArray doesn't contain all the filenames of dataFolder.")
  }
  
  ###CRMAv2 : normalization of CEL files. The same for normal-tumor and single array study
  if (is.null(tags)) {
    ds <- aroma.affymetrix::doCRMAv2(dataFolder, chipType=paste0(chipType),verbose=-1,combineAlleles=FALSE)  
  } else {
    ds <- aroma.affymetrix::doCRMAv2(dataFolder, chipType=paste0(chipType,",",tags),verbose=-1,combineAlleles=FALSE)
  }
  
  
  ###genotype calls: genotypage of the NORMAL data only, not the TUMOR 
  if(genotypeCallsMethod=="naive")
    naiveGenotypeCalls(dataFolder,normalTumorArray,singleStudy,savePlot) 
  
  
  ###tumorBoost: normalization of fraction allele B tumor signal. normal and tumor fracB are required. 
  if(!singleStudy)
    tumorboost(dataFolder,normalTumorArray,savePlot)
 
  return(invisible(NULL)) 
}

#'
#' normalization process for estimating raw copy-numbers and allele B fraction.
#' 
#' @title Normalization process
#' @param dataSetName Name of the data set. If you use architecture=FALSE, the name must correspond to a name of folder in the rawData folder.
#' @param chipType Type of the used chip (e.g. "GenomeWideSNP_6"). If architecture=FALSE, the files of the chip must be contained in the annotationData folder,
#'  if TRUE, they have to be in the "chipTypePath" folder.
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files.
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param dataSetPath (only if createArchitecture=TRUE) Path to the folder containing the CEL files of the data-set.
#' @param chipFilesPath (only if createArchitecture=TRUE) Path to the folder containing all the annotations files for the specified chip type.
#' @param path (only if createArchitecture=TRUE) Path where the architecture should be created (default=".").
#' @param createArchitecture if TRUE, the aroma architecture will be automatically created (default=TRUE). 
#' CEL files of the data and chip files will be copied (not moved).
#' @param savePlot if TRUE, graphics of the CN signal and allele B fraction signal will be saved in the figures/signal folder.
#' @param tags Common tag which appears in the different file names (cdf, ugp, ufl) of the chip. For no tag, use tags=NULL (default = NULL). See details for more information.
#' 
#' @details
#' The following architecture must be used:
#'   <working directory>
#'     +- annotationData/
#'     |  +- chipTypes/
#'     |     +- <chipType>/ <-- must match exactly the name of the CDF file (fullname minus tags)
#'     |        +- CDF file(s) and other annotation (possibly subdirectories)
#'     |
#'     +- rawData/
#'        +- <nameOfDataSet>/
#'           +- <chipType>/ <-- must match exactly a chip type folder under annotationData/
#'              +- CEL files
#'      
#'
#' If you use createArchitecture=TRUE, this function creates this architecture for you and copy your files in the right folders.
#' 
#' The functions will create other folders which contain figures, results of normalization.
#' 
#' If you already have the required architecture, you just have to add your data in the rawData folder with respect to the architecture.
#' 
#' All the cdf chip file names must follow the following rule : <chipType>,<Tags>.cdf
#' 
#' Multiples tags must be separated by  a comma. If there is no tag, the pattern is <chipType>.cdf
#' 
#' @author Quentin Grimonprez
#' 
#' @export
signalPreProcess=function(dataSetName, chipType, normalTumorArray, dataSetPath, chipFilesPath=dataSetPath, path=".", createArchitecture=TRUE, savePlot=TRUE, tags=NULL)
{

  allpkg=TRUE
  if(!require("aroma.affymetrix", quietly=TRUE) )
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
    
  if(!require("aroma.cn", quietly=TRUE) )
  {
    cat("Package not found: aroma.cn. For download it:\n")
    cat("hbInstall(\"aroma.cn\")\n") 
    allpkg=FALSE
  }
  
  if(!allpkg)
    return(NULL)
  
  if(createArchitecture==TRUE)
  {
    actualPath=getwd()
    createArchitecture(dataSetName,chipType,dataSetPath,chipFilesPath,path,TRUE,tags)
    
    #move to the path of the created Architecture
    setwd(path)
  }

  SignalNormalization(dataSetName,chipType,normalTumorArray,"naive",savePlot,tags)
 
  return(invisible(NULL)) 
}


    
  
