#'
#' Extract allele B fraction signals from aroma files. It requires to have executed the normalization process suggested by aroma packages, by using 
#' \link{signalPreProcess} for example.
#'
#' @title Extract allele B fraction signal from aroma files
#' @param dataSetName The name of the data-set folder (it must correspond to a folder name in rawData folder.)
#' @param chromosome A vector containing the chromosomes for which the allele B fraction signal must be extract. 
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param listOfFiles  A vector containing the names of the files in dataSetName folder for which the allele B fraction profiles will be extracted (default is all the files).
#' @param verbose If TRUE print some information (default=TRUE).
#' 
#' @return a list of length the number of chromosomes containing a list of two elements (normal and tumor) containing a data.frame with columns:
#' \describe{
#'   \item{chromosome}{Chromosome of the signal.}
#'   \item{position}{Positions associated with the allele B fraction.}
#'   \item{fracB}{Allele B fraction profiles of selected files; the name of each column is the name of the associated data file name.}
#'   \item{featureNames}{Names of the probes.}
#' }
#' 
#' @details The aroma architecture must be respected. The working directory must contain rawData folder and totalAndFracBData folder.
#' To easily access the names of the files available in a dataset, one can use the \link{getListOfFiles} function.
#'  
#' 
#' @examples 
#' #DO NOT EXECUTE before reading the vignette
#' #fracB=getFracBSignal("data1",5,normalTumorArray)
#' #fracB=getFracBSignal("data2",5)
#'
#' @author Quentin Grimonprez
#'
#' @export 
getFracBSignal=function(dataSetName,chromosome,normalTumorArray,listOfFiles=NULL,verbose=TRUE)
{
  allpkg=TRUE
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

  require(aroma.core)
  require(R.filesets)
  require(R.methodsS3)
  
  if(!("totalAndFracBData"%in%list.files()))
    stop("There is no \"totalAndFracBData\", check if you are in the good working directory or if you have run the signalPreProcess function before.")
  
  ################### check 
  #chromosome
  if(missing(chromosome))
    stop("chromosome is missing.")
  if(!is.numeric(chromosome))
    stop("chromosome must contain integers between 1 and 25.")
  if(sum(unlist(lapply(chromosome,is.wholenumber)))!=length(chromosome))
    stop("chromosome must contain integers between 1 and 25.")
  chromosome=unique(chromosome)
  if(length(chromosome)>25 || length(chromosome)==0)
    stop("chromosome contains too many elements.")
  chromosome=sort(chromosome)
  for(i in 1:length(chromosome))
    if( (chromosome[i]<1) || (chromosome[i]>25) )
      stop("chromosome must contain integers between 1 and 25.")
  
  #dataSetName
  if(missing(dataSetName))
    stop("dataSetName is missing.")
  if(!is.character(dataSetName))
    stop("dataSetName must be the name of an existing folder in rawData.")
  
  #check if we are in a normal-tumor study or in a single array study
  singleStudy=TRUE
  if(missing(normalTumorArray))
  {
    if(verbose)
      cat("No normalTumorArray specified.\n The allele B fraction signal will be extracted for all the specified data.\n")
      #stop("No normalTumorArray specified.\n Youd need to specify a normalTumorArray to extract allele B fraction")
  }
  else
  {
    if(verbose)
     cat("The allele B fraction signal will be extracted for normal and tumor signals. The normalized tumorboost allele B fraction signal will be extracted for tumor signal.\n")
    singleStudy=FALSE
  }
  
  
  ###import dataset to check listOfiles and normalTumorArray
  
  #path where find the fracB data for normal case
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  #for the fracB tumor after tumorboost
  dataSetTumorBoost=paste0(dataSet,",TBN,NGC")
  
  #load fracB
  ds <- aroma.core::AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
  
  
  #check listOfFiles
  pos=c()
  if(is.null(listOfFiles) || missing(listOfFiles))
  {
    listOfFiles=R.filesets::getNames(ds)
    pos=1:length(ds)
  }
  else
  {
    #check the format
    if(!is.vector(listOfFiles) || !is.character(listOfFiles))
      stop("listOfFiles must be a vector of string.")
    listOfFiles=unique(listOfFiles)
    
    #check if all the files of listOfFiles are in the folder
    pos=sapply(listOfFiles,match,R.filesets::getNames(ds))#position of the files of listOfFiles in the folder
    if(length(which(pos>0))!=length(pos))
      stop("Wrong name of files in listOfFiles")
  } 
  
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
    
    #check is the file contains all the file
#     isArrayComplete=sapply(R.filesets::getNames(ds),FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
#     if(sum(isArrayComplete)!=length(isArrayComplete))
#       stop("normalTumorArray doesn't contain all the filenames of dataSetName.")
    
    isArrayComplete=sapply(listOfFiles,FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
    if(sum(isArrayComplete)!=length(isArrayComplete))
      stop("normalTumorArray doesn't contain all the filenames you specified in listOfFiles parameter.")
  }
  
  #if paired study, we keep the name of normal files
  normalFiles=NULL
   if(!singleStudy)
  {
    #if normal-tumor study, we need the tumor and normal files
    
    #we obtain the complementary files
    compFiles=getComplementaryFile(listOfFiles,normalTumorArray)
    allFiles=unique(c(listOfFiles,compFiles))
    
    #get the status ("normal" or "tumor") of each files
    status=getStatus(allFiles,normalTumorArray)
    
    #keep only the normal files
    normalFiles=allFiles[which(status=="normal")]
    
    rm(compFiles,allFiles)
   }   
  
  ########### END CHECK ARGUMENT
  
  ####################
  
  #get names and psoition of the probes
  ugp <- aroma.core::getAromaUgpFile(ds);
  unf <- aroma.core::getUnitNamesFile(ugp);
  
  #get the prefix of SNP probes
  platform <- aroma.core::getPlatform(ugp);
  if (platform == "Affymetrix") 
  {
    require("aroma.affymetrix") || R.methodsS3::throw("Package not loaded: aroma.affymetrix");
    snpPattern <- "^SNP|^S-";
  } 
  else if (platform == "Illumina") 
  {
    snpPattern <- "^rs[0-9]";
  }
  else 
  {
    R.methodsS3::throw("Unknown platform: ", platform);
  }
  
  allFracB=list()
  for(chr in chromosome)
  {
    units <- aroma.core::getUnitsOnChromosome(ugp, chromosome=chr);
  
    unitNames <- aroma.core::getUnitNames(unf,units=units);##names of the probes

      
    #keep the SNP units
    units=units[grep(snpPattern,unitNames)]
    unitNames=unitNames[grep(snpPattern,unitNames)]
  
    posChr <- aroma.core::getPositions(ugp, units=units);#positions of the probes
    #sort signal by position
    indSort=order(posChr)
    
    units=units[indSort]
    unitNames=unitNames[indSort]
    posChr=posChr[indSort]
    
    #get the fracB signal
     if(singleStudy)
     {
       fracB=getFracBSignalSingleStudy(ds,units,pos)
       fracB$tumor=fracB$tumor
     }
     else
     {
      fracB=getFracBSignalPairedStudy(ds,units,normalTumorArray,normalFiles)
      fracB$normal=fracB$normal  
      fracB$tumor=fracB$tumor
      allFracB[[paste0("chr",chr)]]$normal=data.frame(chromosome=rep(chr,length(posChr)),position=posChr,fracB$normal,featureNames=unitNames,check.names=FALSE)
      
     }
    
    allFracB[[paste0("chr",chr)]]$tumor=data.frame(chromosome=rep(chr,length(posChr)),position=posChr,fracB$tumor,featureNames=unitNames,check.names=FALSE)
  }
  
  return(allFracB)
}

###################################################################################################################

#
# @title Get the allele B fraction signal in case of single study
# @param dsC object return by AromaUnitFracBCnBinarySet$byName function
# @param units position to keep
# @param indexOfFiles index of the files to extract
# 
# @return A matrix. each column contains the allele B fraction signal for a different profile.
# 
# @author Quentin Grimonprez
#
getFracBSignalSingleStudy=function(ds,units,indexOfFiles)
{    
  require(aroma.affymetrix)
  require(aroma.cn)
  require(R.filesets)
  
  #reduce to the files from indexOfFiles
  ds=extract(ds,indexOfFiles)
  
  #extract fracB for normal signal
  fracBnormal <- R.filesets::extractMatrix(ds, units=units);
  sampleNames=R.filesets::getNames(ds)
  colnames(fracBnormal)=sampleNames
  
  return(list(tumor=fracBnormal,sampleNames=sampleNames))
}

###################################################################################################################

#
# @title Get fracB signal in case of normal-tumor study
# @param ds object return by AromaUnitFracBCnBinarySet$byName function
# @param units position to keep
# @param normalTumorArray only if you have normal and tumor profile in your data folder. A csv file or a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param normalFiles vector containing the names of normal files to extract
# 
# @return A matrix. each column contains the allele B fraction signal for a different profile.
# 
# @author Quentin Grimonprez
#
getFracBSignalPairedStudy=function(ds,units,normalTumorArray,normalFiles)
{  
  require(aroma.affymetrix)
  require(aroma.cn)
  require(aroma.core)
  require(R.filesets)
  
  #id of normal files
  normalId=sapply(normalFiles,FUN=function(x,names){which(names==x)},R.filesets::getNames(ds))
    
  # Extract selected normal files
  dsFracBnormal <- extract(ds, normalId);  
  
  #extract fracB for normal signal
  fracBnormal <- R.filesets::extractMatrix(dsFracBnormal, units=units);
  
  #folder for tumorboost
  dataSet2=paste0(R.filesets::getFullName(ds),",TBN,NGC")
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dstumor<- aroma.core::AromaUnitFracBCnBinarySet$byName(dataSet2, chipType="*", paths=rootPath);
  
  #names of tumor files
  tumorFiles=getComplementaryFile(normalFiles,normalTumorArray)
  
  #id of tumor files
  tumorId=sapply(tumorFiles,FUN=function(x,names){which(names==x)},R.filesets::getNames(dstumor))

  # Extract selected tumor files
  dsFracBtumor <- extract(dstumor, tumorId);  
  
  #extract fracB for tumor signal
  fracBtumor <- R.filesets::extractMatrix(dsFracBtumor, units=units);
  
  sampleNamesN=R.filesets::getNames(dsFracBnormal)
  sampleNamesT=R.filesets::getNames(dsFracBtumor)

  colnames(fracBnormal)=sampleNamesN
  colnames(fracBtumor)=sampleNamesT
    
  return(list(normal=fracBnormal,tumor=fracBtumor,sampleNamesN=sampleNamesN,sampleNames=sampleNamesT))
}

###################################################################################################################

#'
#' The allele B fraction signal is the ratio between the signal from the allele B and the total signal.
#' The symmetrization of the fraction allele B signal x is : 2*abs(x-0.5).
#'
#' @title symmetrize an allele B fraction signal
#'
#' @param fracB a vector containing an allele B fraction signal.
#'
#' @return a vector containing the symmetrized signal.
#'
#' @examples
#' signalA=abs(rnorm(100))
#' signalB=abs(rnorm(100))
#' signalFracB=signalA/(signalA+signalB)
#' 
#' symFracB=symmetrizeFracB(signalFracB)
#'
#' @export
#'
#' @author Quentin Grimonprez
symmetrizeFracB=function(fracB)
{
  if(missing(fracB))
    stop("fracB is missing.")
  if(!is.numeric(fracB))
    stop("fracB must be a vector of real.")
  #symmetrize the frac B signal: now near 0 it's heterozygous, near 1 homozygous
  return(2*abs(fracB-0.5))
}
