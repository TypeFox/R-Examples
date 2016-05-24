#'
#' Extract genotype calls from aroma files. It requires to have executed the normalization process suggested by aroma packages, by using 
#' \link{signalPreProcess} for example.
#'
#' @title Extract genotype calls from aroma files
#' @param dataSetName The name of the data-set folder (it must correspond to a folder name in rawData folder.)
#' @param chromosome A vector containing the chromosomes for which the genotype call will be extracted. 
#' @param listOfFiles A vector containing the names of the files in dataSetName folder for which the genotype signal will be extracted (default is all the files).
#' @param verbose If TRUE print some information (default=TRUE)
#'
#' @return a list of length the number of chromosomes containing a data.frame with columns:
#' \describe{
#'   \item{chromosome}{Chromosome of the signal.}
#'   \item{position}{Positions associated with the genotype.}
#'   \item{genotype}{Genotype calls corresponding to selected files; the name of each column is the name of the associated data file name.}
#'   \item{featureNames}{Names of the probes.}
#' }
#' 
#' @details The aroma architecture must be respected. The working directory must contain rawData folder and totalAndFracBData folder.
#' To easily access the names of the files available in a dataset, one can use the \link{getListOfFiles} function.
#' 
#' 
#' @examples 
#' #DO NOT EXECUTE before reading the vignette
#' #fracB=getGenotypeCalls("data1",5)
#' 
#' @author Quentin Grimonprez
#' 
#' @export 
getGenotypeCalls=function(dataSetName,chromosome,listOfFiles=NULL,verbose=TRUE)
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
  
  if(!("callData"%in%list.files()))
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
#   singleStudy=TRUE
#   if(missing(normalTumorArray))
#   {
#     if(verbose)
#       cat("No normalTumorArray specified.\n The allele B fraction signal will be extracted for all the specified data.\n")    
#   }
#   else
#   {
#     if(verbose)
#       cat("The allele B fraction signal will be extracted for normal and tumor signals. The normalized tumorboost allele B fraction signal will be extracted for tumor signal.")
#     singleStudy=FALSE
#   }
  
  
  ###import dataset to check listOfiles and normalTumorArray
  
  #path where find the genotype calls data
  rootPath <- "callData";
  rootPath <- Arguments$getReadablePath(rootPath);
  
  genotypeTag <- "NGC";
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  gsN <- aroma.core::AromaUnitGenotypeCallSet$byName(dataSet, tags=genotypeTag, chipType="*");
  
  #check listOfFiles
  pos=c()
  if(is.null(listOfFiles) || missing(listOfFiles))
  {
    listOfFiles=R.filesets::getNames(gsN)
    pos=1:length(gsN)
  }
  else
  {
    #check the format
    if(!is.vector(listOfFiles) || !is.character(listOfFiles))
      stop("listOfFiles must be a vector of string.")
    listOfFiles=unique(listOfFiles)
    
    #check if all the files of listOfFiles are in the folder
    pos=sapply(listOfFiles,match,R.filesets::getNames(gsN))#position of the files of listOfFiles in the folder

    if(length(which(pos>0))!=length(pos))
    {
      phrase="Wrong name of files in listOfFiles.\n Available names are : "
      for(i in R.filesets::getNames(gsN))
        phrase=paste0(phrase,"\n",i)

      stop(phrase)      
    }
  } 

  ########### END CHECK ARGUMENT
  
  ####################
  
  #get names and psoition of the probes
  ugp <- aroma.core::getAromaUgpFile(gsN);
  unf <- aroma.core::getUnitNamesFile(ugp);
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
  else {
    R.methodsS3::throw("Unknown platform: ", platform);
  }
  
  gsN=extract(gsN,pos)
  
  allGen=list()
  for(chr in chromosome)
  {
    units <- aroma.core::getUnitsOnChromosome(ugp, chromosome=chr);
    #get the prefix of SNP probes
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
    
    calls=matrix(nrow=length(posChr),ncol=length(listOfFiles))
    colnames(calls)=R.filesets::getNames(gsN)

    for(i in 1:length(listOfFiles))
    {
      gsNtemp=R.filesets::getFile(gsN,i)
      calls[,i]=gsNtemp[units,1,drop=TRUE]
    }
    calls[calls==0]="BB"
    calls[calls==1]="AB"
    calls[calls==2]="AA"
    
    
    allGen[[paste0("chr",chr)]]=data.frame(chromosome=rep(chr,length(posChr)),position=posChr,calls,featureNames=unitNames)
  }

  return(allGen)
}




