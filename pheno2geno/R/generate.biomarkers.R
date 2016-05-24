#
# generate.biomarkers.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified April, 2013
# first written Mar, 2011
# Contains: generate.biomarkers, pull.biomarkers, selectTopMarker.internal
#           scoreMarker.internal, convertfindBiomarkers.internal
#           splitPheno.internal, selectMarkersUsingMap.internal, 
#           filterGenotypes.internal, filterRow.internal, splitRowSubEM.internal
#

# generate.biomarkers
#
# DESCRIPTION:
#  Function that chooses from the matrix only appropriate markers with specified rules
# PARAMETERS:
#   - population - An object of class population.
#   - threshold - If  pval for gene is lower that this value, we assume it is being diff. expressed.
#   - overlapInd - Number of individuals that are allowed in the overlap.
#   - proportion - Proportion of individuals expected to carrying a certain genotype. 
#   - pProb - Posterior probability threshold used to assugn genotypes.
#   - env - Vector contatining information about environment for each of the individuals in the set (numeric value).
#   - verbose - Be verbose.
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#  An object of class cross
#
generate.biomarkers <- function(population, threshold=0.05, overlapInd = 10, proportion = c(50,50), margin = 15, pProb=0.8, n.cluster=1, env, verbose=FALSE, debugMode=0){
  
  ### checks
  ### check population
  if(missing(population)) stop("Population object not found.\n")
  check.population(population)

  ### check threshold
  if(!is.numeric(threshold))                        stop("threshold should be a numeric value\n")
  if(threshold < 0)                                 stop("threshold should be a positivite value\n")

  ### check overlapInd
  if(!is.numeric(overlapInd))                       stop("overlapInd should be a numeric value\n")
  if(overlapInd < 0)                                stop("overlapInd should be a positivite value\n")

  ### check proportion
  if(any(!is.numeric(proportion)))                  stop("overlapInd should be a numeric value\n")
  if(any(proportion < 1) || sum(proportion) != 100) stop("Proportion should be > 0 and < 100 and sum up to 100\n")

  ### check margin
  if(!is.numeric(margin))                           stop("margin should be a numeric value\n")
  if(margin < 0 || margin > 100)                    stop("margin should be a number between 0 and 100\n")
  
  ### check pProb
  if(!is.numeric(pProb))                            stop("pProb should be a numeric value\n")
  if(pProb < 0 || pProb > 1)                        stop("pProb should be a number between 0 and 1\n")

  if(verbose && debugMode==1) cat("generate.biomarkers starting withour errors in checkpoint.\n")
  s  <- proc.time()
  
  #*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
  s1 <- proc.time()
  if(!is.null(population$annots)){
    population <- generate.biomarkers.internal(population, threshold=threshold, overlapInd=overlapInd, proportion=proportion, 
                                               margin=margin, n.cluster=n.cluster, pProb=pProb, verbose=verbose, debugMode=debugMode)
  }else{
    population <- generate.biomarkers.internal(population, threshold=threshold, overlapInd=overlapInd, proportion=proportion, 
                                               margin=margin, n.cluster=n.cluster, pProb=pProb, verbose=verbose, debugMode=debugMode)
  }
  e1 <- proc.time()
  if(verbose && debugMode==2) cat("Converting phenotypes to genotypes done in:",(e1-s1)[3],"seconds.\n")
  
  #*******RETURNING CROSS OBJECT*******
  e<-proc.time()
  if(verbose) cat("generate.biomarkers done in",(e-s)[3],"seconds.\n")
  if(verbose) cat("Selected",nrow(population$offspring$genotypes$simulated),"markers\n")
  invisible(population)
}

############################################################################################################
#                  *** pull.biomarkers ***
#
# DESCRIPTION:
#  function returning all biomarkers or top marker matching given pattern
# 
# PARAMETERS:
#   population - an object of class population
#   pattern - vector of 0s and 1s (or 0,1,2s)
#   verbose - be verbose
# 
# OUTPUT:
#  vector/matrix
#
############################################################################################################
pull.biomarkers <- function(population,pattern,verbose=FALSE){
  if(missing(population)) stop("No population object provided.\n")
  if(is.null(population$offspring$genotypes$simulated)) stop("Population object doesn't contain de novo genotypes, run findBiomarkers.\n")
  markers <- population$offspring$genotypes$simulated
  if(verbose) cat("Selected",nrow(markers),"markers.\n")
  if(!missing(pattern)){
    if(length(pattern)!=ncol(markers)) stop("Wrong length of the pattern: ",length(pattern)," instead of: ",ncol(markers)," \n")
    if(verbose) cat("Selecting marker best matching given pattern.\n")
    markers <- selectTopMarker.internal(markers,pattern,verbose)
  }
  invisible(markers)
}

############################################################################################################
#                  *** selectTopMarker.internal  ***
#
# DESCRIPTION:
#  function returning all biomarkers or top marker matching given pattern
# 
# PARAMETERS:
#   population - an object of class population
#   pattern - vector of 0s and 1s (or 0,1,2s)
#   verbose - be verbose
# 
# OUTPUT:
#  vector/matrix
#
############################################################################################################
selectTopMarker.internal <- function(markers,pattern,verbose){
  markerPoints <- apply(markers,1,function(x){sum(x==pattern)})
  topMarker    <- rownames(markers)[which.max(markerPoints)]
  if(verbose) cat("Markers best matching pattern:",topMarker,"with identity:",max(markerPoints)/ncol(markers)*100,"%\n")
  invisible(markers[topMarker,])
}

############################################################################################################
#                  *** convertfindBiomarkers.internal ***
#
# DESCRIPTION:
#  function splitting differentially expressed markers into two genotypes
# 
# PARAMETERS:
#   population - object of class population, must contain founders phenotypic data.
#   orderUsing- which map should be used to order markers (default - none)
#     - map_genetic - genetic map
#    - map_physical - physical map
#   treshold - if Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
#   overlapInd - number of individuals that are allowed in the overlap
#   proportion - proportion of individuals expected to carrying a certain genotype 
#   margin - proportion is allowed to varry between this margin (2 sided)
#   verbose - be verbose
#   debugMode - 1: Print our checks, 2: print additional time information 
# 
# OUTPUT:
#  object of class population
#
############################################################################################################
generate.biomarkers.internal <- function(population, threshold, overlapInd, proportion, margin, n.cluster, pProb=0.8, verbose=FALSE, debugMode=0){
  ### initialization
  populationType          <- class(population)[2]
  if(verbose && debugMode==1) cat("generate.biomarkers.internal starting.\n")
  output                  <- NULL
  markerNames             <- NULL 
  
  ### selection step
  if(is.character(population$offspring$phenotypes)){
    offspringFile <- population$offspring$phenotypes
    population$offspring$phenotypes <- NULL
    ### in HT mode markers are read from file line by line and processed on the fly
    population                           <- applyFunctionToFile(offspringFile,sep="\t", header=TRUE, FUN=selectByLine, 
    population=population, threshold=threshold, overlapInd=overlapInd, proportion=proportion, margin=margin, pProb=pProb, n.cluster=n.cluster, verbose=verbose)
    colnames(population$offspring$genotypes$simulated) <- colnames(population$offspring$phenotypes)
    return(population)
  }else{
    ### checking if any of the phenotypes is down/up-regulated
    upRegulatedPhenos       <- selectPhenotypes(population, threshold, 1)
    downRegulatedPhenos     <- selectPhenotypes(population, threshold, 2)
  }
  
  ### removing probes that were selected as both up and down regulated - obsolete as rankprod is not really used any more and t.test will never do that
  if(!is.null(rownames(upRegulatedPhenos)) && !is.null(rownames(downRegulatedPhenos))){
    inupndown               <- which(rownames(upRegulatedPhenos) %in% rownames(downRegulatedPhenos))
    if(length(inupndown)>0) upRegulatedPhenos  <- upRegulatedPhenos[-inupndown,]
  }

  ### if any of the phenotypes is up-regulated - process them
  if(!(is.null(dim(upRegulatedPhenos)))&&(nrow(upRegulatedPhenos)!=0)){
    if(verbose) cat("Selected ",nrow(upRegulatedPhenos),"upregulated potential markers.\n")
    cur                     <- splitPheno.internal(upRegulatedPhenos, overlapInd=overlapInd, proportion=proportion, margin=margin, 
                               pProb=pProb, populationType=populationType, n.cluster=n.cluster, up=TRUE, done=0, left=0, verbose=verbose)
    output                  <- rbind(output,cur)
    
  }else{
    if(verbose) cat("Selected none upregulated potential markers.\n")
  }

  ### if any of the phenotypes is down-regulated - process them
  if(!(is.null(dim(downRegulatedPhenos)))&&(nrow(downRegulatedPhenos)!=0)){
    if(verbose) cat("Selected ",nrow(downRegulatedPhenos),"downregulated potential markers.\n")
    cur                     <- splitPheno.internal(downRegulatedPhenos, overlapInd=overlapInd, proportion=proportion, margin=margin, 
                               pProb=pProb, populationType=populationType, n.cluster=n.cluster, up=FALSE, done=0, left=0, verbose=verbose)

    output                  <- rbind(output,cur)
  }else{
    if(verbose) cat("Selected none downregulated potential markers.\n")
  }
  
  if(verbose) cat("Generated ",nrow(output),"markers.\n")
  ### putting results inside population object
  if(is.null(dim(output))) stop("No markers selected.")

  population$offspring$genotypes$simulated           <- output
  colnames(population$offspring$genotypes$simulated) <- colnames(upRegulatedPhenos)
  invisible(population)
}

### select phenotypes that are suitable for EM algorithm
selectPhenotypes <- function(population, threshold, RPcolumn){
  notNullPhenotypes   <- which(population$founders$RP$pval[,RPcolumn] > 0)        # rank product gives a score for 0 sometimes -> this is below the threshold but these phenotypes wshould not be selected
  belowTreshold       <- which(population$founders$RP$pval[,RPcolumn] < threshold) # phenos diff expressed with pval lower than threshold
  selected            <- belowTreshold[which(belowTreshold%in%notNullPhenotypes)]
  selectedParental    <- population$founders$phenotypes[selected,]
  rownamesOfSelected  <- rownames(selectedParental)
  if(any(rownamesOfSelected == "")) rownamesOfSelected <- rownamesOfSelected[-which(rownamesOfSelected == "")]
  selectedRils        <- population$offspring$phenotypes[rownamesOfSelected,]
  invisible(selectedRils)
}

selectByLine <- function(dataMatrix, population, lineNR, threshold, overlapInd, proportion, margin, pProb, n.cluster){
  result    <- lapply(1:nrow(dataMatrix),selectByLineApply, dataMatrix, population, lineNR, threshold, overlapInd, proportion, margin, pProb, n.cluster)
  resultM   <- do.call(rbind,result)
  population$offspring$genotypes$simulated <- rbind(population$offspring$genotypes$simulated,resultM)
  if(is.null(population$offspring$phenotypes)){
      population$offspring$phenotypes <- rbind(population$offspring$phenotypes,dataMatrix)
  }else if(nrow(population$offspring$phenotypes)<100000){
      population$offspring$phenotypes <- rbind(population$offspring$phenotypes,dataMatrix)
   }
  invisible(population)
}


### select phenotypes that are suitable for EM algorithm in a line by line fashion
selectByLineApply <- function(dataRowNr, dataMatrix, population, lineNR, threshold, overlapInd, proportion, margin, pProb, n.cluster){
  #if(verbose && debugMode==1) cat("selectByLine starting.\n")
  dataRow         <- matrix(dataMatrix[dataRowNr,],1,ncol(dataMatrix))
  dataRowName     <- rownames(dataMatrix)[dataRowNr]
  populationType  <- class(population)[2]
  ### if there is an annotation for that probe - lets use it, if not - do nothing
  if(!is.null(population$annots)){
    ### if there is an annotation data, the name of the probe should be a number corresponding
    # to a number of a row storing information about the probe in annotation matrix
    if(!is.numeric(dataRowName)){
      dataRowName <- as.numeric(dataRowName)
      if(!is.finite(dataRowName)) stop(dataRowName," is not a correct name for a probe!")
    } ### if it is numeric, we still need to check wheter there is a row like that,
    ### if that is true -> select the names of the probe stored there
    if(dataRowName > 0 && dataRowName < nrow(population$annots)) dataRowName <- population$annots[dataRowName,1]
  }
  
  ### is there parental information
  if(!is.null(population$founders$RP$pval)){
    ### is there parental information for a certain probe
    if(dataRowName%in%rownames(population$founders$RP$pval)){
      ### is there parental information for a certain probe
      if(!is.null(population$founders$RP$pval[dataRowName,])){
        ### if it is there but does not pass the threshold - return NULL
        if(!any(population$founders$RP$pval[dataRowName,]>0 & population$founders$RP$pval[dataRowName,]<threshold)) return(NULL)
      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
  }else{
    ### analyse variance of the probe - is it even worth touching by EM
    if(!analyseLineVariance(dataRow,threshold)) return(NULL)
  }
  
  ### split the probe and select [[1]], [[2]] -> info about EM that we cannot store in HT mode
  if((population$founders$RP$pval[dataRowName,2]==0)|(population$founders$RP$pval[dataRowName,1]<population$founders$RP$pval[dataRowName,2])){
    up <-  0
  }else{
    up <- 1
  }
  result        <- splitPheno.internal(dataRow, overlapInd=overlapInd, proportion=proportion, margin=margin, 
                pProb=pProb, populationType=populationType, n.cluster=n.cluster, up=up, done=0, left=0)

  ### if the probe is selected (so result != NULL) return both genotype and phenotype
  if(!is.null(result)){
    ### reformatting as a matrix for easier handling
    result                                         <- matrix(result,1,ncol(result))
    rownames(result)                               <- rownames(dataMatrix)[dataRowNr]
  }
  invisible(result)
}

analyseLineVariance <- function(dataRow,threshold){
  if(any(!is.numeric(dataRow))) invisible(FALSE)
  if(any(is.na(dataRow)))       dataRow <- dataRow[-which(is.na(dataRow))]
  ### code duplication !!! remember to remove it
  half         <- floor(length(dataRow)/2)
  end          <- length(dataRow)
  dataRow      <- sort(dataRow)
  meansToTest  <- c( mean(dataRow[1:half],na.rm=TRUE), 
                     mean(dataRow[2:(half+1)],na.rm=TRUE),
                     mean(dataRow[3:(half+2)],na.rm=TRUE),
                     mean(dataRow[(half+1):end],na.rm=TRUE),
                     mean(dataRow[(half):(end-1)],na.rm=TRUE),
                     mean(dataRow[(half-1):(end-2)],na.rm=TRUE))
  ### end of duplication
  
  ### is the variance passing the threshold?
  res      <- t.test(meansToTest[1:3], meansToTest[4:6])
  if(res$p.val < threshold) return(TRUE)
  invisible(FALSE)
}

# selectPhenotypes
#
# DESCRIPTION:
#  Function that selects offsprings fulfilling the criteria
# PARAMETERS:
#   - population - An object of class population.
#   - threshold - If  pval for gene is lower that this value, we assume it is being diff. expressed.
#
# OUTPUT:
#  A matrix with selected phenotypes
#

mergeEnv.internal <- function(population, genoMatrix){
  ### check if there is anything to merge and if so -> merge
  done    <- NULL
  newGeno <- NULL
  for(probenr in 1:nrow(genoMatrix)){
    probe     <- genoMatrix[probenr,]
    probeNr   <- probe[1]                      # first element of the probe is its number
    probe     <- probe[-1]
    probeID   <- population$annots[probe[1],2] # rows must be matching between annotations and genotypes
    probeName <- population$annots[probe[1],1]
    if(!(probeID %in% done)){                     # if the probe was not yet merged -> we should analyse it
      done        <- c(done,probeID)              # so that we know we have processed it already
      probeNrs    <- which(population$annots[,2]==probeID)
      #cat(probes,":",length(probes),"\n")
      
      if(length(probeNrs)>1){                       # we have more than one probe with the same ID -> merging
        cat("Mergining:",length(probeNrs),"probes with ID:",probeID,"\n")
        probes        <- genoMatrix[probeNrs,]
        # for each of the positions we set a consensus genotype
        consensusGeno <- round(apply(probes,2,mean,na.rm=TRUE))
      }
    }
    newGeno <- rbind(newGeno,c(probeName,probe))
  }
  invisible(newGeno)
}

############################################################################################################
#                  *** splitPheno.internal ***
#
# DESCRIPTION:
#  subfunction of convertfindBiomarkers.internal, splitting children markers using founders mean values
# 
# PARAMETERS:
#   offspring - matrix of up/down regulated genes in offspring
#   founders - matrix of up/down regulated genes in parents
#   overlapInd - Number of individuals that are allowed in the overlap
#   proportion - Proportion of individuals expected to carrying a certain genotype 
#   margin - Proportion is allowed to varry between this margin (2 sided)
#   groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
#   up - 1 - genes up 0 - down regulated
# 
# OUTPUT:
#  list containg genotype matrix and names of selected markers
#
############################################################################################################
splitPheno.internal <- function(offspring, overlapInd, proportion, margin, populationType, up, n.cluster=1, pProb=0.8, done=0, left=0, verbose=FALSE){
  output           <- NULL
  s                <- proc.time()
  printedProc      <- NULL
  #cl               <- makeCluster(getOption("cl.cores", n.cluster))
  results          <- lapply(1:nrow(offspring), splitPheno.Apply, offspring=offspring, overlapInd=overlapInd, proportion=proportion, margin=margin, pProb=pProb, up=up, populationType = populationType, verbose=verbose)
  #stopCluster(cl)
  output           <- do.call(rbind,results)
  invisible(output)
}

splitPheno.Apply <- function(x, offspring, overlapInd, proportion, margin, pProb, up, populationType, verbose){
  cur <- splitPhenoRowEM.internal(offspring[x,], overlapInd, proportion, margin, pProb, up, populationType, verbose)
  if(!(is.null(cur))){
    output           <- matrix(cur, 1, length(cur))
    rownames(output) <- rownames(offspring)[x]
    return(output)
  }
  return(NULL)
}

############################################################################################################
#                  *** splitPhenoRowEM.internal ***
#
# DESCRIPTION:
#  subfunction of splitRow.internal, splitting one row using EM algorithm
# 
# PARAMETERS:
#   x - name of currently processed row
#   offspring - matrix of up/down regulated genes in offspring
#   founders - matrix of up/down regulated genes in parents
#   overlapInd - Number of individuals that are allowed in the overlap
#   proportion - Proportion of individuals expected to carrying a certain genotype 
#   margin - Proportion is allowed to varry between this margin (2 sided)
#   groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
#   up - 1 - genes up 0 - down regulated
# 
# OUTPUT:
#  genotype row
#
############################################################################################################
splitPhenoRowEM.internal <- function(x, overlapInd, proportion, margin, pProb=0.8, up=1, populationType, verbose=FALSE){
  nrDistributions <- length(proportion)
  result          <- rep(NA,length(x))
  EM              <- NULL
  idx             <- which(!(is.na(x)))
  idw             <- length(which((is.na(x))))
  y               <- x[idx]
  if(populationType == "f2"){
    minimalObsRequired <- ceiling(0.67*length(x)) # three normals - 2/3 of data is not NA
    if(length(y) < minimalObsRequired) stop("Too little observations for accurate EM fitting! At least: ",minimalObsRequired," observations required!")
  }else{
    minimalObsRequired <- ceiling(0.5*length(x)) # two normals - 1/2 of data is not NA
    if(length(y) < minimalObsRequired) stop("Too little observations for accurate EM fitting! At least: ",minimalObsRequired," observations required!")
  }
  ### EM
  tryCatch({
    aa <- tempfile()
    sink(aa)
    EM <- normalmixEM(y, k=nrDistributions, lambda= proportion, maxrestarts=1, maxit = 300, fast=FALSE)
  },
  error= function(err){
    warning(paste("ERROR in splitPhenoRowEM.internal while running EM:  ",err))
    sink()            # sink if errored -> otherwise everything is sinked into aa file
    # file is not removed -> contains output that may help with debugging
  },
  finally={
    sink()
    file.remove(aa) # no error -> close sink and remove unneeded file
  })
  if(filterRow.internal(EM$lambda,proportion,margin)){
    if(populationType == "f2"){
      genotypes <- c(1:5)
      if(up==0) genotypes <- c(3,2,1,5,4)
      for(i in (1:length(y))){
        if(any(EM$posterior[i,]>pProb)){
          result[idx[i]] <- genotypes[which.max(EM$posterior[i,])]
        }else if((EM$posterior[i,1]+EM$posterior[i,2])>pProb){
          result[idx[i]] <- genotypes[4]
        }else if((EM$posterior[i,2]+EM$posterior[i,3])>pProb){
          result[idx[i]] <- genotypes[5]
        }else{
          result[idx[i]] <- NA
        }
      }
    }else{
      genotypes <- c(1:2)
      if(up==0) genotypes <- c(2,1)
      for(i in (1:length(y))){
        result[idx[i]] <- NA
        if(any(EM$posterior[i,]>pProb)) result[idx[i]] <- genotypes[which.max(EM$posterior[i,])]
      }
    }
  }else{
    result <- NULL
  }
  if(!is.null(result)) if((sum(is.na(result))-idw)>overlapInd) result <- NULL
  invisible(result)
}


############################################################################################################
#                  *** filterRowSub.internal ***
#
# DESCRIPTION:
#   subfunction of filterGenotypes.internal, filtering one row
# 
# PARAMETERS:
#   genotypeRow - currently processed row
#   overlapInd - Number of individuals that are allowed in the overlap
#   proportion - Proportion of individuals expected to carrying a certain genotype 
#   margin - Proportion is allowed to varry between this margin (2 sided)
# 
# OUTPUT:
#  boolean
#
############################################################################################################
filterRow.internal<- function(lambda, proportion, margin){
  if(length(lambda) != length(proportion)) return(FALSE)
  for(i in 1:length(lambda)){
    if((lambda[i]>((proportion[i]+margin/2)/100)) || (lambda[i]<((proportion[i]-margin/2)/100))) return(FALSE)    #TODO use abs and merge the 2 similar IF conditions
  }
  return(TRUE)
}

