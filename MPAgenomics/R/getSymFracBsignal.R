#'
#' Extract symmetrized allele B fraction signals from aroma files. It requires to have executed the normalization process suggested by aroma packages, by using 
#' \link{signalPreProcess} for example.
#'
#' @title Extract symmetrized allele B fraction signal from aroma files
#' @param dataSetName The name of the data-set folder (it must correpond to a folder name in rawData folder.)
#' @param chromosome A vector with the chromosomes for which the symetrized signal will be extracted. 
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param file The name of the file in dataSetName to extract.
#' @param verbose If TRUE, print some informations.
#' 
#' @return a list of length the number of chromosome containing a data.frame with columns:
#' \describe{
#'   \item{chromosome}{chromosome corresponding to the signal.}
#'   \item{position}{Positions associated to the allele B fraction.}
#'   \item{fracB}{One column named by the data file name. It contains the symmetrized allele B fraction signal for the specified profile.}
#'   \item{featureNames}{Names of the probes.}
#' }
#' 
#' @details The aroma architecture must be respected. The working directory must contain rawData folder and totalAndFracBData folder.
#' To easily access the names of the files available in a dataset, one can use the \link{getListOfFiles} function.
#' 
#' @examples 
#' #DO NOT EXECUTE
#' #fracB=getSymFracBSignal("data1",5,normalTumorArray)
#' #fracB=getSymFracBSignal("data2",5)
#'
#' @author Quentin Grimonprez
#'
#' @export 
getSymFracBSignal=function(dataSetName,file,chromosome,normalTumorArray,verbose=TRUE)
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
  #singleStudy=TRUE
  if(missing(normalTumorArray))
  {
    if(verbose)
      #cat("No normalTumorArray specified.\n The allele B fraction signal will be extracted for all the specified data.\n")    
      stop("No normalTumorArray specified.\n Youd need to specify a normalTumorArray to extract and symmetrize allele B fraction")
  }
#   else
#   {
#     if(verbose)
#       cat("The allele B fraction signal will be extracted for normal and tumor signals. The normalized tumorboost allele B fraction signal will be extracted for tumor signal.")
#     singleStudy=FALSE
#   }
  
  
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
  if(missing(file))
    stop("file is missing.")
  else
  {
    #check the format
    if(!is.character(file))
      stop("file must be a string.")
    if(length(character)!=1)
      stop("file must be a string.")
      
    #check if all the files of listOfFiles are in the folder
    pos=match(file,R.filesets::getNames(ds))#position of the files of listOfFiles in the folder
    if(is.na(pos))
      stop("Wrong name of file.")

  }

  
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
    
    #check is the file contains all the file
    #     isArrayComplete=sapply(R.filesets::getNames(ds),FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
    #     if(sum(isArrayComplete)!=length(isArrayComplete))
    #       stop("normalTumorArray doesn't contain all the filenames of dataSetName.")
    
    if(!(file%in%c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor))))
      stop("normalTumorArray doesn't contain all the filenames you specified in listOfFiles parameter.")
#   }
  
  #if paired study, we keep the name of normal files
  normalFiles=NULL
#   if(!singleStudy)
#   {
    #if normal-tumor study, we need the tumor and normal files
    
    #we obtain the complementary files
    compFiles=getComplementaryFile(file,normalTumorArray)
    allFiles=unique(c(file,compFiles))
    
    #get the status ("normal" or "tumor") of each files
    status=getStatus(allFiles,normalTumorArray)
    
    #keep only the normal files
    normalFiles=allFiles[which(status=="normal")]

    rm(compFiles,allFiles)

#   }   


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
    R.methodsS3::throw("Unknown platform: ", platform);
  

  symFracB=list()
  for(chr in chromosome)
  {
    units <- aroma.core::getUnitsOnChromosome(ugp, chromosome=chr);
    unitNames <- aroma.core::getUnitNames(unf,units=units);##names of the probes
    
    
    #keep the SNP units
    units=units[grep(snpPattern,unitNames)]
    unitNames=unitNames[grep(snpPattern,unitNames)]
    
    posChr <- aroma.core::getPositions(ugp, units=units);#positions of the probes
    #sort signal by position
    indSort=sort(posChr,index.return=TRUE)$ix
    
    posChr=posChr[indSort]
    units=units[indSort]
    unitNames=unitNames[indSort]
    
    #get the genotype calls for 1 chr
    geno=getGenotypeCalls(dataSetName,chromosome=chr,listOfFiles=normalFiles,verbose=FALSE)
    geno=geno[[paste0("chr",chr)]]
    
    ind=which(geno[,3]=="AB")
    
    posChr=posChr[ind]
    units=units[ind]
    unitNames=unitNames[ind]
    
    fracB=getFracBSignalPairedStudy(ds,units,normalTumorArray,normalFiles)
    fracB$normal=symmetrizeFracB(fracB$normal)
    fracB$tumor=symmetrizeFracB(fracB$tumor)
    symFracB[[paste0("chr",chr)]]$normal=data.frame(chromosome=rep(chr,length(posChr)),position=posChr,fracB$normal,featureNames=unitNames)
    symFracB[[paste0("chr",chr)]]$tumor=data.frame(chromosome=rep(chr,length(posChr)),position=posChr,fracB$tumor,featureNames=unitNames)
  }

  return(symFracB)
}


