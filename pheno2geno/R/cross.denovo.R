#
# cross.denovo.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Mar, 2011
# Contains: orderChromosomes, majorityRule.internal, mergeChromosomes.internal
#           switchChromosomes.internal, removeChromosomes.internal 
#           removeChromosomesSub.internal
#

# cross.denovo
#
# DESCRIPTION:
#  Ordering chromosomes using genetic/physical map and majority rule
# OUTPUT:
#  An object of class cross
#
cross.denovo <- function(population, n.chr, map=c("none","genetic","physical"), comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation,majorityOfMarkers), 
assignFunction=c(assignMaximumNoConflicts,assignMaximum), reOrder=TRUE, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0, ...){
  #checks
  if(missing(population)) stop("Provide a population object.\n")
  check.population(population)

  if(missing(n.chr))     stop("Provide number of expected chromosomes.\n")
  if(!is.numeric(n.chr)) stop("Number of expected chromosomes must be a numeric value.\n")
  if(n.chr < 1)          stop("Number of expected chromosomes must be a positive value.\n")
  if("noParents" %in% population$flag) n.chr <- n.chr*2 

  map              <- match.arg(map)
  comparisonMethod <- defaultCheck.internal(comparisonMethod,"comparisonMethod",4,sumMajorityCorrelation)
  assignFunction   <- defaultCheck.internal(assignFunction,"assignFunction",2,assignMaximumNoConflicts)

  cross            <- cross.denovo.internal(population,n.chr,verbose=TRUE,debugMode=2)

  if(length(cross$geno)<=1){
    cat("Selected cross object contains too little chromosomes to assign them, returning it.")
    return(cross)
  }
  if(verbose && debugMode==1) cat("cross.denovo starting withour errors in checkpoints.\n")

  if(map=="none"){
    if(reOrder){
      cross <- formLinkageGroups(cross, reorgMarkers=TRUE, ...)
      cross <- reduceChromosomesNumber(cross, n.chr)
      if("noParents" %in% population$flags) cross <- mergeInverted(cross, class(population)[2])
      if(use.orderMarkers) cross <- orderMarkers(cross, use.ripple=TRUE, verbose=TRUE)
      return(cross)
    }else{
      assignment        <- names(cross$geno)
      names(assignment) <- names(cross$geno)
      return(assignment)
    }
  }
  s1 <- proc.time()
  if(map=="genetic"){
    if(is.null(population$maps$genetic==NULL)) stop("No genetic map provided in population$maps$genetic\n")
    originalMap <- population$maps$genetic
  }else if(map=="physical"){
    if(is.null(population$maps$physical)) stop("No physical map provided in population$maps$physical\n")
    originalMap <- population$maps$physical
  }
  cross <- formLinkageGroups(cross,reorgMarkers=TRUE,...)
  cross <- reduceChromosomesNumber(cross, n.chr)
  if("noParents" %in% population$flags){
    cross <- mergeInverted(cross, class(population)[2])
  }
  if(use.orderMarkers){
      cross <- orderMarkers(cross,use.ripple=TRUE,verbose=TRUE)
  }
  chromToChromArray <- comparisonMethod(cross, originalMap, population)
  e1 <- proc.time()
  if(verbose){cat("Calculating correlation matrix done in:",(e1-s1)[3],"seconds.\n")}
  assignment <- assignFunction(chromToChromArray)
  if(!reOrder){
    if(verbose)cat("Returning new ordering vector.\n")
    invisible(assignment)
  }else{
    ordering <- assignChrToMarkers(assignment,cross)
    if(verbose)cat("Applying new ordering to the cross object.\n")
    if(use.orderMarkers){
        if(verbose)cat("Ordering markers inside the cross object\n")
        s0 <- proc.time()
        nmarkersPerChr <- nmar(cross)
        nChr <- length(nmarkersPerChr)
        cross <- reorganizeMarkersWithin(cross,ordering)
        for(i in 1:nChr){
          cross <- orderMarkers(cross, use.ripple=TRUE, chr=i, verbose=TRUE)
          e1 <- proc.time()
          if(i < nChr){
            te <- ((e1-s0)[3]/sum(nmarkersPerChr[1:i]))*sum(nmarkersPerChr[(i+1):nChr])
          }else{ 
            te <- 0 
          }
          if(verbose) cat("Done ordering chromosome",i,"/",nChr,"Time remaining:",te,"seconds.\n")
        }
        e0 <- proc.time()
        if(verbose && debugMode==2)cat("Ordering markers inside the cross object done in:",(e0-s0)[3],"seconds.\n")
    }else{
      cross <- reorganizeMarkersWithin(cross,ordering)
    }
  cross <- replace.map(cross,est.map(cross))
  print(nchr(cross))
  invisible(cross)
  }
}


############################################################################################################
#                  *** assignMaximum ***
#
# DESCRIPTION:
#  function returning for rows or cols - which element is having max value
# OUTPUT:
#  vector of numerics
############################################################################################################
assignMaximum <- function(x, use = 2){ apply(x,use,which.max) }

############################################################################################################
#                  *** assignMaximumNoConflicts ***
#
# DESCRIPTION:
#  function return
# OUTPUT:
#  vector of numerics
############################################################################################################
assignMaximumNoConflicts <- function(x, use = 2){
  assignment <- assignMaximum(x,use)
  notYetAssigned <- as.numeric(names(assignment)[which(!(names(assignment)%in%assignment))])
  while(any(duplicated(assignment))){
    duplicated_ones <- assignment[(duplicated(assignment))]
    for(duplication in duplicated_ones){
      need_to_decide <- which(assignment == duplication)
      best_fitting <- names(which.max(apply(x,use,max)[need_to_decide]))
      need_to_decide <- need_to_decide[-which(names(need_to_decide)==best_fitting)]
      #### we need to work on the next line, probably some while loop here, because if you call function
      #### few times(I mean all the functions on the way not only this one, it's getting better and better!:P
      assignment[as.numeric(names(need_to_decide))] <- notYetAssigned[1]
      notYetAssigned <- as.numeric(names(assignment)[which(!(names(assignment)%in%assignment))])
      }
    }
  invisible(assignment)
}

mergeInverted <- function(cross,populationType){
  chr.correlations           <- matrix(0,nchr(cross),nchr(cross))
  mar.correlations           <- cor(pull.geno(cross),use="pair")
  rownames(mar.correlations) <- markernames(cross)
  colnames(mar.correlations) <- markernames(cross)
  colnames(chr.correlations) <- 1:nchr(cross)
  rownames(chr.correlations) <- 1:nchr(cross)
  for(chr in 1:nchr(cross)){
    markers1 <- colnames(cross$geno[[chr]]$data)
    for(chr2 in 1:nchr(cross)){
      markers2 <- colnames(cross$geno[[chr2]]$data)
      chr.correlations[chr,chr2] <- mean(mar.correlations[markers1,markers2],na.rm=TRUE)
    }
  }
  corPairs <- t(apply(chr.correlations,1,function(x){names(sort(x))}))
  colnames(corPairs) <- 1:nchr(cross)
  done <- NULL
  new <- 1:nchr(cross)
  cur <- 1
  ordering <- vector(mode="numeric",length=sum(nmar(cross)))
  names(ordering) <- markernames(cross)
  for(chr in 1:nchr(cross)){
    if(!(chr %in% done)){
      done <- c(done,chr)
      chr2 <- selectChr(chr,corPairs,done)
      if(!is.null(chr2)){
        done <- c(done,chr2)
        MinVal <- chr.correlations[chr,chr2]
        if(abs(MinVal)>0){
            markers <- colnames(cross$geno[[chr]]$data)
            print(markers)
            markers <- c(markers,colnames(cross$geno[[chr2]]$data))
            ordering[markers] <- new[cur]
            cur <- cur + 1
            if(MinVal<0){
              RefVal <- chr.correlations[chr,1]
              if(RefVal < 0){
                cross <- invertChromosome.internal(cross,chr,populationType)
              }else{
                cross <- invertChromosome.internal(cross,chr2,populationType)
              }
            }
        }else{
          if(chr.correlations[chr,1] < chr.correlations[chr2,1]){
            cross <- removeChromosomes(cross,chr)
          }else{
            cross <- removeChromosomes(cross,chr)
          }
        }
      }
    }
  }
  cross <- reorganizeMarkersWithin(cross, ordering)
  invisible(cross)
}

selectChr <- function(chr,corPairs,done){
  for(i in colnames(corPairs)){
    if(!(corPairs[chr,i] %in% done)){
      return(corPairs[chr,i])
    }
  }
  return(NULL)
}

invertChromosome.internal <- function(cross,chr,populationType){
  cat("--- Inverting chromosome:",chr,"---\n")
  if(populationType=="f2"){
    for(i in chr){
      if(i %in% 1:nchr(cross)){
        geno_ <- cross$geno[[i]]$data
        if(any(geno_==5)) cross$geno[[i]]$data[which(geno_==5)] <- 4
        if(any(geno_==4)) cross$geno[[i]]$data[which(geno_==4)] <- 5
        if(any(geno_==3)) cross$geno[[i]]$data[which(geno_==3)] <- 1
        if(any(geno_==1)) cross$geno[[i]]$data[which(geno_==1)] <- 3
      }else{
        stop("Chromosome: ",i," is not in the set!\n")
      }
    }
  }else{
    for(i in chr){
      if(i %in% 1:nchr(cross)){
        geno_ <- cross$geno[[i]]$data
        if(any(geno_==2)) cross$geno[[i]]$data[which(geno_==2)] <- 1
        if(any(geno_==1)) cross$geno[[i]]$data[which(geno_==1)] <- 2
      }else{
        stop("Chromosome: ",i," is not in the set!\n")
      }
    }
  }
  invisible(cross)
}

############################################################################################################
#                                           *** majorityRule ***
#
# DESCRIPTION:
#   Subfunction of assignChromosomes, for every chromosome in cross for every marker checks the marker it is
#   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
#   new chromosomee, old chromosome with most markers with high correlation is assigned. 
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
majorityCorrelation <- function(cross, originalMap, population, verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)
  chrMaxCorMarkers <- t(apply(abs(genotypesCorelationMatrix),2,function(r){c(originalMap[rownames(genotypesCorelationMatrix)[which.max(r)],1],max(r))}))
  ordering <- NULL
  chromToChromArray <- matrix(0,length(unique(originalMap[,1])),nrOfChromosomesInCross)
  rownames(chromToChromArray) <- unique(originalMap[,1])
  colnames(chromToChromArray) <- 1:nrOfChromosomesInCross
  for(i in 1:nrOfChromosomesInCross){
    markersFromCurChrom <-colnames(cross$geno[[i]]$data)
    #### USING ONLY MARKERS WITH COR HIGHER THAN 0.5
    markersCorWithCurChrom <- chrMaxCorMarkers[markersFromCurChrom,]
    if(length(markersFromCurChrom)>1){
        #markersHighlyCorWithCurChrom <- markersCorWithCurChrom[which(abs(markersCorWithCurChrom[,2])>0.5),]
        chromosomesCorWithCurrent <- table(markersCorWithCurChrom[,1])
        #bestCorChrom <- as.numeric(names(chromosomesCorWithCurrent)[which.max(chromosomesCorWithCurrent)])
        for(j in names(chromosomesCorWithCurrent)){
          chromToChromArray[j,i] <- chromosomesCorWithCurrent[j]
        }
    }else if(length(markersFromCurChrom)==1){
        #markersHighlyCorWithCurChrom <- markersCorWithCurChrom[which(abs(markersCorWithCurChrom[,2])>0.5),]
        chromToChromArray[markersCorWithCurChrom[1],i] <- markersCorWithCurChrom[2]
    }else{
         bestCorChrom <- NA
    }
    
  }
  invisible(chromToChromArray)
}



############################################################################################################
#                                          ** sumMajorityRule***
#
# DESCRIPTION:
#   Subfunction of assignChromosomes, for every chromosome in cross for every marker checks the marker it is
#   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
#   new chromosomes one or more of chromosomes from old map will be represented. Function sums correlations for
#   each pair of those and for every new chromosomes assigns old chromosome with highest cumulative cor.
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
sumMajorityCorrelation <- function(cross, originalMap, population, verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)
  chrMaxCorMarkers <- t(apply(abs(genotypesCorelationMatrix),2,function(r){c(originalMap[rownames(genotypesCorelationMatrix)[which.max(r)],1],max(r))}))
  chromToChromArray <- matrix(0,length(unique(originalMap[,1])),nrOfChromosomesInCross)
  rownames(chromToChromArray) <- unique(originalMap[,1])
  colnames(chromToChromArray) <- 1:nrOfChromosomesInCross
  for(i in 1:nrOfChromosomesInCross){
    markersFromCurChrom <-colnames(cross$geno[[i]]$data)
    markersHighlyCorWithCurChrom <- chrMaxCorMarkers[markersFromCurChrom,]
    if(length(markersFromCurChrom)>1){
        correlatedChrom <- as.numeric(names(table(markersHighlyCorWithCurChrom[,1])))
        for(j in correlatedChrom){
              currentSum <- sum(markersHighlyCorWithCurChrom[which(markersHighlyCorWithCurChrom[,1]==j),2])
              chromToChromArray[j,i] <- currentSum
        }
    }else{
       chromToChromArray[markersHighlyCorWithCurChrom[1],i] <- markersHighlyCorWithCurChrom[2]
    }
  }
  invisible(chromToChromArray)
}

############################################################################################################
#                                         *** correlationRule ***
#
# DESCRIPTION:
#  Subfunction of assignChromosomes, assigning chromosome from new map to old ones using mean corelation\
#  between their markers.
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
meanCorrelation <- function(cross,originalMap,population,verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)
  chromToChromArray <- matrix(0,length(unique(originalMap[,1])),nrOfChromosomesInCross)
  rownames(chromToChromArray) <- unique(originalMap[,1])
  colnames(chromToChromArray) <- 1:nrOfChromosomesInCross
  for(i in 1:nrOfChromosomesInCross){
    markersFromCurNewChrom <-colnames(cross$geno[[i]]$data)
    for(j in unique(originalMap[,1])){
        markersFromCurOldChrom <- rownames(originalMap)[which(originalMap[,1]==j)]
        currentCorMatrix <- abs(genotypesCorelationMatrix[markersFromCurOldChrom,markersFromCurNewChrom])
        chromToChromArray[j,i] <- sum(apply(currentCorMatrix,1,mean))
    }
  }
  return(chromToChromArray)
}

############################################################################################################
#                                         *** correlationRule ***
#
# DESCRIPTION:
#  Subfunction of assignChromosomes, assigning chromosome from new map to old ones using mean corelation\
#  between their markers.
# 
# PARAMETERS:
#  cross - object of class cross
#  originalMap - map from population object
#  genotypesCorelationMatrix - gene correlation matrix (from map2mapCorrelationMatrix function)
#  verbose - be verbose
# 
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
majorityOfMarkers <- function(cross,originalMap,population,verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  chromToChromArray <- matrix(0,length(unique(originalMap[,1])),nrOfChromosomesInCross)
  rownames(chromToChromArray) <- unique(originalMap[,1])
  colnames(chromToChromArray) <- 1:nrOfChromosomesInCross
  for(i in 1:nrOfChromosomesInCross){
    markersFromCurNewChrom <-colnames(cross$geno[[i]]$data)
    for(j in unique(originalMap[,1])){
        markersFromCurOldChrom <- rownames(originalMap)[which(originalMap[,1]==j)]
        chromToChromArray[j,i] <- sum(markersFromCurOldChrom%in%markersFromCurNewChrom)
    }
  }
  return(chromToChromArray)
}


############################################################################################################
#                  *** cross.denovo.internal ***
#
# DESCRIPTION:
#   function to create new map and save it in cross object
# 
# PARAMETERS:
#   population - object of class population
#   orde - object of class population
#   n.chr - expected number of linkage groups
#   use - expected number of linkage groups
#  verbose - be verbose
#
# OUTPUT:
#  an object of class cross
#
#
############################################################################################################
cross.denovo.internal<- function(population,  n.chr,  use=c("rf","geno"), verbose=FALSE, debugMode=0){
  if(missing(n.chr)) stop("n.chr in an obligatory parameter")
  if(missing(population)) stop("no population object provided")
  use <- match.arg(use)
  check.population(population)
  if(is.null(population$offspring$genotypes$simulated)){
    stop("no simulated genotypes in population object, first use generate.biomarkers!\n")
  }
  #*******SAVING CROSS OBJECT*******
  s1 <- proc.time()
  ### creation of the cross
  tryCatch({
    aa <- tempfile()
    sink(aa)
    cross <- genotypesToCross.internal(population,"simulated",verbose=verbose,debugMode=debugMode)
  },
  error= function(err){
    stop(paste("ERROR in cross.saturate while creating cross:  ",err))
    sink()            # sink if errored -> otherwise everything is sinked into aa file
    # file is not removed -> contains output that may help with debugging
  },
  finally={
    sink()
    file.remove(aa) # no error -> close sink and remove unneeded file
  })
  e1 <- proc.time()
  if(verbose && debugMode==2)cat("saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
  
  #*******CREATING NEW MAP*******
  s1 <- proc.time()
  cross <- assignLinkageGroups(cross,n.chr,use="rf")
  e1 <- proc.time()
  if(verbose && debugMode==2)cat("New map created in:",(e1-s1)[3],"seconds.\n")
  
  invisible(cross)
}
