#
# readFiles.r
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Nov, 2012
# first written Nov, 2011
# Contains: read.population mapMarkers.internal, correctRowLoc.internal
#           probesLocation.internal, orrectRowGff.internal
#

#  read.population
#
# DESCRIPTION:
#  Reads geno/phenotypic files into R environment into special object.
# PARAMETERS:
#   - offspring - Core used to specify names of children phenotypic ("offspring_phenotypes.txt") and genotypic ("offspring_genotypes.txt") files.
#   - founders - Core used to specify names of founders phenotypic ("founders_phenotypes.txt") file.
#   - map - Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files.
#   - founders_groups - specify founders groups
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#   An object of class population 
#
read.population <- function(offspring = "offspring", founders = "founders", map = "map", foundersGroups, populationType = c("riself", "f2", "bc", "risib"),
  readMode = c("normal","HT"), threshold=0.05, verbose = FALSE, debugMode = 0, n.cluster=1, ...){

  ### checks
  populationType <- match.arg(populationType)
  readMode <- match.arg(readMode)
  
  
  ### file names
  fileFoundersPheno  <- paste(founders,"_phenotypes.txt",sep="")
  fileOffspringPheno <- paste(offspring,"_phenotypes.txt",sep="")
  fileOffspringGeno  <- paste(offspring,"_genotypes.txt",sep="")
  fileAnnotations    <- paste(offspring,"_annotations.txt",sep="")
  fileMapPhys        <- paste(map,"_physical.txt",sep="")
  fileMapGen         <- paste(map,"_genetic.txt",sep="")

  ### initializing  
  s <- proc.time()
  if(verbose && debugMode==1) cat("read.population starting.\n")
  population <- NULL
  
  ### offspring phenotypic file
  if(!file.exists(fileOffspringPheno)){
    stop("No phenotype file for offspring: ",fileOffspringPheno," this file is essential, you have to provide it\n")
  }else{
    if(verbose)  cat("File:",fileOffspringPheno,"found and will be processed.\n")
    if(readMode == "normal"){
      population                      <- applyFunctionToFile(fileOffspringPheno, population, sep="\t", header=TRUE, FUN=normalModeReading)
      if(verbose) cat("   + Read in",nrow(population$offspring$phenotypes),"phenotypes\n")
    }else{
      population$offspring$phenotypes <- fileOffspringPheno
    }
  }

  ### founders phenotypic file
  if(!file.exists(fileFoundersPheno)){
    ### simulate data if there is no file
    if(verbose)cat("No phenotype file for founders: ",fileFoundersPheno,". Founder phenotypes will be simulated.\n")
    
    if(readMode == "normal"){
      population                         <- simulateParentalPhenotypes(population, population$offspring$phenotypes, populationType)
    }
    ### if the mode is HT, we don't simulate founders but just judge the variance in the offspring while converting phenotypes to genotypes
    
  }else{
    ### read the file if present
    if(verbose)  cat("File:",fileFoundersPheno,"found and will be processed.\n")
    
    ### check if there is an information about the founders groups
    if(missing(foundersGroups)) stop("No information about founders groups provided.\n")
    
    ### founders groups should be a sequence of 0s and 1s
    if(any(foundersGroups!=0 && foundersGroups!=1)) stop("Founders groups attribute is incorrect.\n")
    population                        <- applyFunctionToFile(fileFoundersPheno,population,sep="\t", header=TRUE, verbose=verbose, dataGroups=foundersGroups, threshold=threshold, n.cluster=n.cluster, FUN=tTestByLine)
    population$founders$groups        <- foundersGroups
    if(length(population$founders$groups)   != ncol(population$founders$phenotypes)) stop("Founders groups attribute is incorrect.\n")
    if(verbose) cat("   + Read in",nrow(population$founders$phenotypes),"phenotypes\n")
  }
  
  class(population) <- c("population",populationType)
  
  ### annotations file
  population <- readSingleFile(population, fileAnnotations, "annotations", verbose=verbose, header=TRUE)

  ### offspring genotypic file
  population <- readSingleFile(population, fileOffspringGeno, "offspring$genotypes", verbose=verbose, header=TRUE)
  
  ### physical map
  population <- readSingleFile(population, fileMapPhys, "maps$physical", verbose=verbose, header=FALSE)
  
  ### genetic map
  population <- readSingleFile(population, fileMapGen, "maps$genetic", verbose=verbose, header=FALSE)

  ### TO BE REMOVED when generate.biomarkers is corrected
  population$sliceSize <- 5000
  #**********FINALIZING FUNCTION*************
  e <- proc.time()
  if(verbose && debugMode==2) cat("read.population finished after",(e-s)[3],"seconds.\n")
  invisible(population)
}

#  readSingleFile
#
# DESCRIPTION:
#  Reads single geno/phenotypic file into R environment into special object.
# PARAMETERS:
#   - population - an object of class population
#   - filename - name of the file that will be processed
#   - fileType - type of the file that will be processed (information about what type of data it contains)
# OUTPUT:
#   An object of class population 
#
readSingleFile   <- function(population, filename, fileType, verbose=FALSE, ...){
  if(file.exists(filename)){
    if(verbose)  cat("File:",filename,"found and will be processed.\n")
    dataRead     <- read.table(filename,sep="\t", row.names=1, ...)
    population   <- add.to.population(population, dataRead, fileType)
    doCleanUp.internal()
  }else{
    if(verbose) cat("File:",filename,"not found.\n")
  }
  invisible(population)
}

#  applyFunctionToFile
#
# DESCRIPTION:
#  Reads a file line by line, applying a function to each of the lines.
# PARAMETERS:
#   - filename - name of the file that will be processed
#   - header - does the file contain header
#   - sep - separator of the values in the file
#   - FUN - function to be applied to the lines of the file 
#   - ... - parameters passed to FUN
# OUTPUT:
#   A matrix with values from the file.
#
applyFunctionToFile <- function(filename, population, header=TRUE, sep="\t", chunkSize = 50000, FUN, verbose=FALSE, ...){
  filePointer <- file(filename,"r")
  if(header){
    headerLine <- readLines(filePointer, n=1)
    header     <- unlist(strsplit(headerLine,sep))
  }
  lineNR       <- 0
  ### reading the first non-header line
  curLines <- readLines(filePointer, n=chunkSize)
  while(length(curLines) > 0){
    lineNR                  <- lineNR + chunkSize
    if(verbose && lineNR%%10000==0) cat("      processing line:",lineNR,"\n")
    curLineSplitted                 <- strsplit(curLines,sep)
    curLineSplittedMatrix           <- do.call(rbind,curLineSplitted)
    curRownames                     <- curLineSplittedMatrix[,1]
    curLineSplittedMatrix           <- curLineSplittedMatrix[,-1]
    if(chunkSize==1) curLineSplittedMatrix <- matrix(curLineSplittedMatrix,1,length(curLineSplittedMatrix))
    rownames(curLineSplittedMatrix) <- curRownames
    curLineSplittedMatrix           <- apply(curLineSplittedMatrix,c(1,2),as.numeric)

    ### changing it into a matrix for easier handling
    
    ### if there is header, use it as colnames
    if(!is.null(header)){
      if(length(header)!= ncol(curLineSplittedMatrix)){
        stop("Incorect length of line: ",lineNR," it is: ",ncol(curLineSplittedMatrix)," instead of: ",length(header),"\n")
      }else{
        colnames(curLineSplittedMatrix) <- header
      }
    }
    
    ### execute the function specified by user and rbind results
    population <-  FUN(curLineSplittedMatrix, population, lineNR, ...)
    curLines <- readLines(filePointer, n=chunkSize)
  }
  close(filePointer)
  invisible(population)
}

#  tTestByLine
#
# DESCRIPTION:
#  T.test a single line of a file (formated as a matrix with one row)
# PARAMETERS:
#   - dataMatrix - matrix containing data (with one row)
#   - dataGroups - specifing which columns of the matrix belong to which group
#   - threshold - threshold for pval
# OUTPUT:
#   An input matrix, if the pval is below the threshold. Otherwise NULL
#

tTestByLine <- function(dataMatrix, population, lineNR, dataGroups, threshold, n.cluster=1){
  #cl                             <- makeCluster(getOption("cl.cores", n.cluster))
  #result                         <- parApply(dataMatrix, 1, tTestByLineApply, lineNR, dataGroups, threshold )
  result                         <- apply(dataMatrix, 1, tTestByLineApply, lineNR, dataGroups, threshold )
  #stopCluster(cl)
  resultM                        <- do.call(rbind,result)
  population$founders$phenotypes <- rbind(population$founders$phenotypes,resultM)
  invisible(population)
}

tTestByLineApply <- function(dataRow, lineNR, dataGroups, threshold){
  if(length(dataRow)!=length(dataGroups)) stop("Incorrect line ",lineNR," \n")
  group0     <- as.numeric(dataRow[which(dataGroups == 0)])
  group1     <- as.numeric(dataRow[which(dataGroups == 1)])
  if(length(group0)<3 || length(group1)<3) stop("Not enough observations to perform the t.test.\n")
  res        <- t.test(group0, group1)

  if(res$p.value < threshold) return(dataRow)
  invisible(NULL)
}

#  normalModeReading
#
# DESCRIPTION:
#  Returns the input object without changing it
# PARAMETERS:
#   - dataMatrix - matrix containing data (with one row)
# OUTPUT:
#   An input matrix.
#
normalModeReading <- function(dataMatrix, population, lineNR){
  population$offspring$phenotypes <- checkAndBind(population$offspring$phenotypes,dataMatrix,lineNR)
  invisible(population)
}


checkAndBind <- function(dataMatrix, toBind, lineNR){
 ### if the population object is not empty then we need to check if we can put it into phenotype matrix 
  if(!is.null(dataMatrix)){
    ### can we rbind it?
    if(ncol(toBind) != ncol(dataMatrix)){
      stop("Incorect length of line: ",lineNR," it is: ",ncol(toBind)," instead of: ",ncol(dataMatrix),"\n")
    }
  }### if the object is still empty - we need to fill it
  dataMatrix <- rbind(dataMatrix,toBind)
  invisible(dataMatrix)
}

#  simulateParentalPhenotypes
#
# DESCRIPTION:
#  Simulating founders phenotypes based on offspring data
# PARAMETERS:
#   - population - an object of class population
#   - populationType - breeeding scheme used in the population
#   - offspringPhenotypes - matrix containing phenotypes of the offspring
# OUTPUT:
#   An object of class population
#
simulateParentalPhenotypes <- function(population, offspringPhenotypes, populationType){
  cat("No founders phenotype data provided, it will be simulated!\n")
  half     <- floor(ncol(offspringPhenotypes)/2)
  end      <- ncol(offspringPhenotypes)
  founders <- t(apply(offspringPhenotypes, 1, function(x){
    x      <- sort(x)
    c( mean(x[1:half],na.rm=TRUE), mean(x[2:(half+1)],na.rm=TRUE), mean(x[3:(half+2)],na.rm=TRUE),
       mean(x[(half+1):end],na.rm=TRUE), mean(x[(half):(end-1)],na.rm=TRUE), mean(x[(half-1):(end-2)],na.rm=TRUE))
  }))
  population$flags           <- c(population$flags,"noParents")
  population                 <- add.to.populationSub.internal(population, founders, "founders", populationType=populationType)
  population$founders$groups <- c(0,0,0,1,1,1)
  return(population)
}


############################################################################################################
#                  *** mapMarkers.internal ***
#
# DESCRIPTION:
#  removes from matrix1 cols or rows, which are not present in second (coparing using col/rownames)
# 
# PARAMETERS:
#   expressionMatrix1, expressionMatrix2 - matrices with data of any type
#   mapMode - 1 - map rows, 2 - map cols
#   verbose - Be verbose
#   debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#  object of class population 
#
############################################################################################################
mapMarkers.internal <- function(expressionMatrix1, expressionMatrix2, mapMode=2, verbose=FALSE, debugMode=0){
  if(mapMode==1) {
    nrRows <- nrow(expressionMatrix1)
    ### warnings when names are mismatching
    if(verbose && debugMode==2)if(nrRows!=nrow(expressionMatrix2)){
      cat("Following markers will be removed:\n")
      cat(paste(rownames(expressionMatrix1)[which(!(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)))],"\n"))
    }
    ### mapping itself
    expressionMatrix1 <- expressionMatrix1[which(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)),]
    if(verbose) cat("Because of names mismatch,",nrRows-nrow(expressionMatrix1),"markers were removed, run function with verbose=T debugMode=2 to print their names out.\n")
  }
  else if(mapMode==2){
    nrCols <- ncol(expressionMatrix1)
    ### warnings when names are mismatching
    if(verbose && debugMode==2)if(nrCols!=ncol(expressionMatrix2)){
      cat("Following individuals will be removed:\n")
      paste(colnames(expressionMatrix1)[which(!(colnames(expressionMatrix1) %in% colnames(expressionMatrix2)))],"\n")
    }
    ### mapping itself
    expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
    if(verbose) cat("Because of names mismatch,",nrCols-ncol(expressionMatrix1),"individuals were removed, run function with verbose=T debugMode=2 to print their names out.\n")
  }
  invisible(expressionMatrix1)
}
