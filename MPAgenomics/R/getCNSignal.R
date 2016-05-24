#' 
#' Extract copy-number signals from aroma files. It requires to have executed the normalization process suggested by aroma packages, by using 
#' \link{signalPreProcess} for example.
#'
#' @title Extract copy-number signal from aroma files
#' @param dataSetName The name of the data-set folder (it must correspond to a folder name in rawData folder.).
#' @param chromosome A vector containing the chromosomes for which the signal will be extracted. 
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files.
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param onlySNP If TRUE, only the copy-number for SNPs positions will be returned (default=FALSE).
#' @param listOfFiles A vector containing the names of the files in dataSetName folder for which the copy-number profiles will be extracted (default is all the files).
#' @param verbose If TRUE print some information (default=TRUE).
#' 
#' @return a list of length the number of chromosomes containing a data.frame with columns:
#' \describe{
#'   \item{chromosome}{Chromosome of the signal.}
#'   \item{position}{Positions associated with the copy-number.}
#'   \item{copynumber}{Copy number profiles of selected files; the name of each column is the name of the associated data file name.}
#'   \item{featureNames}{Names of the probes.}
#' }
#'  
#' @details The aroma architecture must be respected. The working directory must contain rawData folder and totalAndFracBData folder.
#' To easily access the names of the files available in a dataset, one can use the \link{getListOfFiles} function.
#' 
#' @examples 
#' #DO NOT EXECUTE before reading the vignette
#' #C=getCopyNumberSignal("data1",5,normalTumorArray,TRUE)
#' #C=getCopyNumberSignal("data2",5,onlySNP=TRUE)
#'
#' @author Quentin Grimonprez
#'
#' @export 
getCopyNumberSignal=function(dataSetName,chromosome,normalTumorArray,onlySNP=FALSE,listOfFiles=NULL,verbose=TRUE)
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
  require(R.methodsS3)
  
  #TODO verifier comportement list.files() sur une autre machine
  list.files()
  if(!("totalAndFracBData"%in%list.files()))
    stop("There is no \"totalAndFracBData\", check if you are in the good working directory or if you have run the signalPreProcess function before.")

  ################### check 
  #onlySNP
  if(!is.logical(onlySNP))
    stop("onlySNP must be either TRUE or FALSE.")
  
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
      cat("No normalTumorArray specified.\n The copy-number signal will be extracted for all the specified data and the median of the dataset will be used as reference.\n")
  }
  else
  {
    if(verbose)
      cat("The copy-number signal will be extracted for all the specified data and the copy-number signal of normal DNA will be used as reference.\n")
    singleStudy=FALSE
  }
  
  #################### import the dataSet to have the name of all the files
  
  #path where find the CN data
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  #load CN
  dsC <- aroma.core::AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
  #print(dsC);
  
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
    
    #check is the normalTumorArray file contains all the files
#     isArrayComplete=sapply(R.filesets::getNames(dsC),FUN=function(name,nameOfFiles){name%in%nameOfFiles},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
#     if(sum(isArrayComplete)!=length(isArrayComplete))
#       stop("normalTumorArray doesn't contain all the filenames of dataSetName.")
    
    if(!missing(normalTumorArray))
    {
      isArrayComplete=sapply(listOfFiles,FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
      if(sum(isArrayComplete)!=length(isArrayComplete))
        stop("normalTumorArray doesn't contain all the filenames you specified in listOfFiles parameter.")
    }
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
  
  
  #get names and position of the probes
  ugp <- aroma.core::getAromaUgpFile(dsC);
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
  
  allCN=list()
  for(chr in chromosome)
  {
    units <- aroma.core::getUnitsOnChromosome(ugp, chromosome=chr);
    unitNames <- aroma.core::getUnitNames(unf,units=units);##names of the probes
    if(onlySNP)
    {
      #keep the SNP units
      units=units[grep(snpPattern,unitNames)]
      unitNames=unitNames[grep(snpPattern,unitNames)]
    }
    posChr <- aroma.core::getPositions(ugp, units=units);#positions of the probes
    #sort signal by position
    indSort=order(posChr)
    
    units=units[indSort]
    unitNames=unitNames[indSort]
    posChr=posChr[indSort]
    
    
    #get the normalized CN signal
    if(singleStudy)
      C=getCopyNumberSignalSingleStudy(dsC,units,chromosome,pos,ugp)
    else
      C=getCopyNumberSignalPairedStudy(dsC,units,normalTumorArray,chromosome,normalFiles,ugp)
    

    allCN[[paste0("chr",chr)]]=data.frame(chromosome=rep(chr,length(posChr)),position=posChr,C$CN[,,drop=FALSE],featureNames=unitNames)
    colnames(allCN[[paste0("chr",chr)]])[-c(1,2,ncol(allCN[[paste0("chr",chr)]]))]=colnames(C$CN)
    
  }

  return(allCN)
}

##########################################################################################################################################

#
# @title Get CN signal in case of single study
# @param dsC object return by AromaUnitTotalCnBinarySet$byName function
# @param units position to keep
# @param chromosome chromosome to extract
# @param indexOfFiles position (in dsC) of the files to extract
# 
# @return A list of 2 elements :
# \describe{
#   \item{CN}{A matrix. Each column contains the CN signal for a different profile.}
#   \item{sampleNames}{Names of the extracted signal.}
# }
# 
# @author Quentin Grimonprez
getCopyNumberSignalSingleStudy=function(dsC,units,chromosome,indexOfFiles,ugp)
{    
  require(aroma.affymetrix)
  require(aroma.core)
  require(R.filesets)
  require(R.oo)
  
  #compute the median
  ceR <- aroma.core::getAverageFile(dsC, verbose=0)
  
  #reduce to the files from indexOfFiles
  dsC=extract(dsC,indexOfFiles)
  
  # Extract total CNs
  C <- R.filesets::extractMatrix(dsC,units=units);
  
  #CN median for units
  thetaR <- R.filesets::extractMatrix(ceR,units=units)
  
  sampleNames=R.filesets::getNames(dsC)
  
  #the ploidy depends of the gender for the cromosome 23 and 24
  if(chromosome==23 || chromosome ==24)
  {
    gender=findGender(R.oo::getName(dsC),indexOfFiles,ugp)
    ploidy=rep(0,length(gender))
    ploidy[gender=="XY"]=1
    if(chromosome==23)
      ploidy[gender=="XX"]=2
    else
      ploidy[gender=="XX"]=NA
    
    C <- sapply(1:ncol(C),FUN=function(x){ploidy[x]*C[,x]/thetaR});
    #BUG FIX : colnames of C are NULL after sapply. reafectation of the proper colnames
    colnames(C)=sampleNames
  }
  else
  {
    #normalization
    C <- apply(C,2,FUN=function(x){2*x/thetaR});
  }

  
  return(list(CN=C,sampleNames=sampleNames))
}

##########################################################################################################################################

#
# @title Get CN signal in case of normal-tumor study
# @param dsC object return by AromaUnitTotalCnBinarySet$byName function
# @param units position to keep
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param chromosome chromosome to extract
# @param normalFiles names of the normal files
# 
# @return A list of 2 elements :
# \describe{
#   \item{CN}{A matrix. Each column contains the CN signal for a different profile.}
#   \item{sampleNames}{Names of the extracted signal.}
# }
# 
# @author Quentin Grimonprez
getCopyNumberSignalPairedStudy=function(dsC,units,normalTumorArray,chromosome,normalFiles,ugp)
{  
  require(aroma.affymetrix)
  require(aroma.core)
  require(R.filesets)
  require(R.oo)
  
  #id of normal files
  normalId=sapply(normalFiles,FUN=function(x,names){which(names==x)},R.filesets::getNames(dsC))
  
  #names of tumor files
  tumorFiles=getComplementaryFile(normalFiles,normalTumorArray)
  
  #id of tumor files
  tumorId=sapply(tumorFiles,FUN=function(x,names){which(names==x)},R.filesets::getNames(dsC))

  # Extract normal
  dsPairCNormal <- extract(dsC, normalId);  
  
  # Extract tumor
  dsPairCTumor <- extract(dsC, tumorId);  
  
  #tumor CN
  C <- R.filesets::extractMatrix(dsPairCTumor, units=units);
  
  #normalCN
  Cnormal <- R.filesets::extractMatrix(dsPairCNormal, units=units);
  sampleNames=R.filesets::getNames(dsC)[tumorId]
  #the ploidy depends of the gender for the cromosome 23 and 24
  if(chromosome==23 || chromosome ==24)
  {
    gender=findGender(R.oo::getName(dsC),normalId,ugp)
    ploidy=rep(0,length(gender))
    ploidy[gender=="XY"]=1#male:chr23=chrX, chr24=chrY, the ploidy is 1
    if(chromosome==23)
      ploidy[gender=="XX"]=2#female:chr23=chrXX, the ploidy is 2
    else
      ploidy[gender=="XX"]=NA#female:chr24=nothing, the ploidy is NA
    
    #normalize by the normal sample
    C <- C/Cnormal
    #apply the ploidy
    C <- sapply(1:ncol(C),FUN=function(x){ploidy[x]*C[,x]});
    #BUG FIX : colnames of C are NULL after sapply. reafectation of the proper colnames
    colnames(C)=sampleNames
  }
  else
  {
    #normalization of CN
    C <- 2*C/Cnormal;
  }

  sampleNames=R.filesets::getNames(dsC)[tumorId]
  
  return(list(CN=C,sampleNames=sampleNames))
}


##########################################################################################################################################

#find gender of a profile CEL file
#
# @param name of the dataset (it must correpond to a folder name in rawData folder.).
# @param indexFileToExtract index of the files to test
#
# @return a vector containing the gender ("XX" or "XY") of the specified files
# @author Quentin Grimonprez
findGender=function(dataSetName,indexFileToExtract,ugp)
{
  require(aroma.cn)
  require(aroma.core)
  require(R.filesets)
  
  #path where find the data
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  #get the allele B fraction
  dsF <- aroma.core::AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
        
  # Identify units on ChrX and ChrY 
  units23 <- aroma.core::getUnitsOnChromosome(ugp, 23);
  units24 <- aroma.core::getUnitsOnChromosome(ugp, 24);
    
  gender=lapply(indexFileToExtract,FUN=function(index)
  {
    df=R.filesets::getFile(dsF,index)
  
    beta23=df[units23,1,drop=TRUE]
    beta24=df[units24,1,drop=TRUE]
    
    genderTemp <- aroma.cn::callXXorXY(beta23, beta24, adjust=1.5, from=0, to=1);
        
  })
  
  return(unlist(gender))
}

##########################################################################################################################################

# check if a number is an integer
# @param x number
# @param tol tolerance
#
# @ return TRUE if the number is an integer, FALSE else
# @author Quentin Grimonprez
is.wholenumber=function(x, tol = .Machine$double.eps^0.5)  
{
  #if(!is.double(x))
  #  return(FALSE)
  
  abs(x - round(x)) < tol
}

##########################################################################################################################################

# get the normal files associated to a tumor files (or vice-versa)
#
# @param listOfFiles vector containing names of iles
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
#
# @return a vector of the same size as listOfFiles containing the complementary files
# @author Quentin Grimonprez
getComplementaryFile=function(listOfFiles,normalTumorArray)
{
  sapply(listOfFiles,FUN=function(names,data)
  {
    ind=match(names,data$normal)
    if(!is.na(ind))
    {
      return(as.character(data$tumor[ind]))
    }
    else
    {
      ind=match(names,data$tumor)
      if(!is.na(ind))
        return(as.character(data$normal[ind]))
      else
        return(NA)
    }
  },normalTumorArray);
}

##########################################################################################################################################

# get the status of a file ("normal" or "tumor")
#
# @param listOfFiles vector containing names of iles
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
#
# @return a vector containing the status of each files of listOfFiles
# @author Quentin Grimonprez
getStatus=function(listOfFiles,normalTumorArray)
{
  sapply(listOfFiles,FUN=function(names,data)
  {
    if(names%in%data$normal)
    {
      return("normal")
    }
    else
    {
      if(names%in%data$tumor)
        return("tumor")
      else
        return(NA)
    }
  },normalTumorArray);
}

##########################################################################################################################################

# get the normal files associated to a tumor files (or vice-versa)
#
# @param listOfFiles vector containing names of iles
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
#
# @return a vector of the same size as listOfFiles containing the complementary files
# @author Quentin Grimonprez
getId=function(listOfFiles,normalTumorArray)
{
  unlist(lapply(listOfFiles,FUN=function(names,data)
  {
    ind=match(names,data$normal)
    if(!is.na(ind))
    {
      return(ind)
    }
    else
    {
      ind=match(names,data$tumor)
      if(!is.na(ind))
        return(ind)
      else
        return(NA)
    }
  },normalTumorArray));
}

# Return a matrix which match normal and tumor files
#
# @param listOfFiles vector containing names of iles
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
#
# @return a matrix of size [number of tumor files, 2]
# @author Samuel BLANCK
getNormalTumorMatrix=function(listOfFiles,normalTumorArray)
{
  status = getStatus(listOfFiles,normalTumorArray)
  
  #id of the tumor files
  indTumor = which(status=="tumor")
  
  #names of complementary normal files
  normFiles=getComplementaryFile(listOfFiles[indTumor],normalTumorArray)
    
  #id of the complementary normal files
  indNormal=sapply(normFiles,FUN=function(x,names){which(names==x)}, listOfFiles)
  
  #matching matrix (1st column=normal files indices, 2nd column=tumor files indices)
  normalTumorMatrix=matrix(ncol=2, nrow=length(indTumor))
  normalTumorMatrix[,1]=indNormal
  normalTumorMatrix[,2]=indTumor
    
  return (normalTumorMatrix)
}