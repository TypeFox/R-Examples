#
# cross.saturate.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified October, 2013
# first written Mar, 2011
# Contains: cross.saturate, rearrangeMarkers, bestCorelated.internal
#           map2mapCorrelationMatrix.internal, map2mapImage
#

# cross.saturate
#
# DESCRIPTION:
#  Saturate an existing genetic map by adding markers derived from expression
# OUTPUT:
#  An object of class cross
#
cross.saturate <- function(population, cross, map=c("genetic","physical"), placeUsing=c("qtl","correlation"), flagged = c("remove","warn","ignore"), threshold=3, chr, env, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0){

  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
  populationType <- class(population)[2]
  
  map            <- match.arg(map)
  placeUsing     <- match.arg(placeUsing)
  flagged        <- match.arg(flagged)
  
  if(missing(env)) env <- rep(1,ncol(population$offspring$phenotypes))
  if(length(env)!=ncol(population$offspring$phenotypes)) stop("Incorrect environmental vector!\n")
  
  if(!is.numeric(threshold)||is.na(threshold)) stop("Please provide correct threshold")
  if(threshold<0) stop("Threshold needs to be > 0")

  if(placeUsing=="correlation" && threshold >= 5) cat("WARNING: threshold too high, few new markers will be selected\n")
  if(placeUsing=="qtl" && threshold >= 20)        cat("WARNING: threshold too high, few new markers will be selected\n")

  ### if the cross object is not provided -> we can create it from original genotypes and map stored in population object
  if(missing(cross)){
    ### checking if we are  able to recreate the original cross object
    if(is.null(population$offspring$genotypes$real))              stop("No original genotypes in population$offspring$genotypes$real, load them in using add.to.population")
    if(is.null(population$offspring$genotypes$simulated))         stop("No genotype data in population$offspring$genotypes$simulated, run generate.biomarkers first")
  }else{
    population   <- set.geno.from.cross(cross,population,map)
    population   <- scan.qtls(population,map,env=env)
  }
  startTime1     <- proc.time()
  
  ### creation of the cross
  tryCatch({
    aa <- tempfile()
    sink(aa)

    populationSubset                      <- population
    populationSubset$offspring$phenotypes <- matrix(0, 5, ncol(population$offspring$phenotypes))
    colnames(populationSubset$offspring$phenotypes) <- colnames(population$offspring$phenotypes)
    rownames(populationSubset$offspring$phenotypes) <- 1:5

    cross  <- genotypesToCross.internal(populationSubset,"simulated",verbose=verbose,debugMode=debugMode)
    cross$pheno <- t(population$offspring$phenotypes)
    file.remove(aa) # no error -> close sink and remove unneeded file
  },
  error   = function(err){ stop(paste("ERROR in cross.saturate while creating cross:  ",err)) },
  finally = { sink() })
  
  markerNames <- rownames(population$offspring$genotypes$simulated)
  lodNames <- rownames(population$offspring$genotypes$qtl$lod)

  if(!(all(markerNames %in% lodNames))) stop("QTL scan results don't match with simulated genotypes, please, run scan.qtls function")
  if(!(all(lodNames %in% markerNames))) stop("QTL scan results don't match with simulated genotypes, please, run scan.qtls function")
  
 if(map=="genetic"){
    population      <- matchMarkers(population, population$maps$genetic, mapType="genetic")
    originalMap     <- population$maps$genetic
  }else{
    population      <- matchMarkers(population, population$maps$physical, mapType="physical")
    originalMap     <- population$maps$physical
  }
  
  nrOfOriginalMarkers  <- nrow(population$offspring$genotypes$real)

  ### saturating only a subset of chromosomes
  if(missing(chr)){
    if(verbose) cat("Saturating all the chromosomes in the set.\n")
    chr           <- unique(originalMap[,1])
  }else{
    availableChr  <- unique(originalMap[,1])
    if(any(!(chr%in%availableChr))) stop("Incorrect chr parameter!\n")
    if(verbose) cat("Saturating chromosomes:\n",paste(chr,",",sep=""),"\n")
  }

  #*******ENRICHING ORIGINAL MAP*******
  cross      <- rearrangeMarkers(cross, population, populationType, originalMap, threshold, placeUsing,
                             flagged, env, addMarkers=TRUE, chr, verbose=verbose, debugMode=debugMode)
  envMarkers <- cross$envMarkers # in case order.markers are used, this info will be erased
  epiMarkers <- cross$epiMarkers
  endTime1   <- proc.time()
  if(verbose && debugMode==2) cat("++ Saturation of the original map done in:",(endTime1-startTime1)[3],"seconds.\n")
  
  #*******ORDERING NEW MAP*******
  if(use.orderMarkers){
    if(verbose) cat("Ordering markers inside the cross object\n")
    startTime1  <- proc.time()
    tryCatch({
      aa <- tempfile()
      sink(aa)
      cross          <- orderMarkers(cross,use.ripple=FALSE,verbose=TRUE)
      file.remove(aa) # no error -> close sink and remove unneeded file
    },
    error= function(err){ stop(paste("ERROR in cross.saturate while ordering markers in cross:  ",err)) },
    finally={ sink() })

    endTime1    <- proc.time()
    if(verbose && debugMode==2)cat("++ Saving data into cross object done in:",(endTime1-startTime1)[3],"seconds.\n")
  }

  nrOfNewMarkers <- sum(nmar(cross))-  nrOfOriginalMarkers 
  percentageOfSaturation   <- (nrOfNewMarkers /nrOfOriginalMarkers )*100
  if(verbose){
    cat("\n=== Saturation statistics:\n")
    cat("Number of original markers:      ", nrOfOriginalMarkers ,"\n")
    cat("Number of inserted markers:      ", nrOfNewMarkers ,"\n")
    cat("Saturation (% of markers added): ", percentageOfSaturation,"\n\n")
  }
  cross$envMarkers <- envMarkers
  cross$epiMarkers <- epiMarkers
  invisible(cross)
}

matchMarkers <- function(population, map, mapType=c("genetic","physical")){
  mapType          <- match.arg(mapType)
  matchingMarkers  <- which(rownames(population$offspring$genotypes$real)%in%rownames(map))
  
  if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
  
  if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
    population$offspring$genotypes$real   <- population$offspring$genotypes$real[matchingMarkers,]
    map                                   <- map[rownames(population$offspring$genotypes$real),]

    cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
  }
  if(mapType=="genetic"){
    population$maps$genetic  <- map
  }else{
    population$maps$physical <- map
  }
  invisible(population)
}

###########################################################################################################
#                                    *** rearrangeMarkers ***
#
# DESCRIPTION:
#   ordering chromosomes using genetic/physical map and corelation rule
# OUTPUT:
#  object of class cross
############################################################################################################
rearrangeMarkers <- function(cross, population, populationType, originalMap, threshold=3, placeUsing, flagged, env, addMarkers=FALSE, chr, verbose=FALSE, debugMode=0){
  ### addMarkers will be removed - obsolete

  if(verbose) cat("Original map contains",max(originalMap[,1]),"chromosomes.\n")
  
  if(placeUsing=="qtl"){
    markersOutput      <- bestQTL.internal(cross,population,threshold,flagged,env,verbose,debugMode)
    markersNewPostions <- markersOutput[[1]]
    markerNames        <- markersOutput[[2]]
    envMarkers         <- markersOutput[[3]]
    epiMarkers         <- markersOutput[[4]]
  }else{
    markersNewPostions <- bestCorelated.internal(cross,population,originalMap,threshold,verbose)
  }
  markersToBeRemoved <- rownames(markersNewPostions) #removing phenotypes used for marker creation from phenotype lists
  ### markersNewPostions - matrix: rows - markers, columns: chr where marker is mapping - position on chr - LOD/cor score
  if(verbose) cat("Selected:\n\t",nrow(markersNewPostions),"markers for further analysis\n\n")

  ### creating a cross object with empty genotypes
  returncross      <- cross
  returncross$geno <- vector(length(unique(originalMap[,1])), mode="list")

  ### removing phenotypes used as markers
  if(any(colnames(returncross$pheno)%in%markersToBeRemoved)){
    returncross$pheno <- returncross$pheno[,-which(colnames(returncross$pheno)%in%markersToBeRemoved)]
  }
  
  ### names of original markers
  originalNames <- rownames(originalMap)[which(rownames(originalMap) %in% rownames(population$offspring$genotypes$real))]
  
  for(x in 1:length(returncross$geno)){
    originalNamesFromCurChr <- rownames(originalMap)[which(originalMap[,1]==x)]
    originalNamesSelected   <- originalNamesFromCurChr[which(originalNamesFromCurChr %in% originalNames)]
    originalPositions       <- originalMap[originalNamesSelected,2]
    if(x %in% chr){
      newnames                  <- rownames(markersNewPostions)[which(markersNewPostions[,1]==x)]
      ### remove new markers that are already in the map
      if(any(newnames%in%originalNamesSelected)) newnames <- newnames[-which(newnames%in%originalNamesSelected)]
      newnamesSelected          <- NULL
      newpositions              <- NULL
      positions                 <- cbind(newnames,markersNewPostions[newnames,2],markersNewPostions[newnames,3])
      ### are there markers mapping to the same position
      for(position in unique(positions[,2])){
        mappingMarkers  <- which(positions[,2]==position)
        newpositions    <- c(newpositions, position)
        if(length(mappingMarkers)==1){
          newnamesSelected  <- c(newnamesSelected,positions[mappingMarkers,1])
        }else{
          ### if more than one marker maps to a certain position - select the one with the highest (QLT/correlation) score
          selM              <- positions[mappingMarkers,1]
          bestM             <- which.max(as.numeric(positions[selM,3]))
          if(length(bestM) > 1) bestM <- bestM[1]
          newnamesSelected  <- c(newnamesSelected,selM[bestM])
        }
      }
    }else{
      newnamesSelected <- NULL
      newpositions <- NULL
    }
    ### constructing the map
    if(verbose){
      cat("Chromosome",x,"\n")
      if(x %in% chr){ cat("\tSelected:",length(newnamesSelected),"new and",length(originalNamesSelected),"original markers.\n")
      }else{ cat("\tSelected:",length(originalNamesSelected),"original markers.\n")}
    }
    returncross$geno[[x]]$data           <- insertMarkers.internal(pull.geno(cross)[,newnamesSelected],newpositions,t(population$offspring$genotypes$real[originalNamesSelected,]), originalPositions, env, populationType)
    newmap                               <- c(as.numeric(newpositions),originalPositions) # new and old positions together
    names(newmap)                        <- c(newnamesSelected,originalNamesSelected)     # new and old names together
    colnames(returncross$geno[[x]]$data) <- c(newnamesSelected,originalNamesSelected)     # new and old names together
    newmap                               <- sort(newmap)                                  # markers on the map must be order based on their position
    returncross$geno[[x]]$data           <- returncross$geno[[x]]$data[,names(newmap)]    # the same order in map and in data matrix
    returncross$geno[[x]]$map            <- c(newmap)
  }
  names(returncross$geno) <- 1:nchr(returncross) ### setting chrnames
  
  ### class of chromosomes set to autosomes - required by r/qtl
  for(i in 1:length(returncross$geno)){
    class(returncross$geno[[i]]) <- "A"
  }
  returncross$envMarkers <- envMarkers
  returncross$epiMarkers <- epiMarkers
  invisible(returncross)
}

###
insertMarkers.internal <- function(newgeno,newpositions,oldgeno,originalPositions,env,populationType){
  if(length(newgeno)<1){ return(oldgeno) }
  toRmv <- NULL    #markers to be removed
  toInv <- NULL    #markers to be inverted
  if(is.null(dim(newgeno))){ newgeno <- as.matrix(newgeno) }  #Does this do anything ? If there is no dim how would as.matrix figure it out then ?
  if(is.null(dim(oldgeno))){ oldgeno <- as.matrix(oldgeno) }  #Does this do anything ? If there is no dim how would as.matrix figure it out then ?

  for(i in 1:length(newpositions)){
    distance <- abs(originalPositions-as.numeric(newpositions[i]))        # calculating distance of new markers to original ones
    curCor   <- cor(newgeno[,i],oldgeno[,which.min(distance)],use="pair") # correlation with the closest original marker
    if(abs(curCor) < 0.1){     ### TODO? let user set this value
      ### if correlation to the closest original markers is that low the markers should be removed
      toRmv <- c(toRmv,i)
    }else if(curCor < (-0.4)){ ### TODO? let user set this value
      ### if correlation to the closest original markers is negative, the marker should be inverted
      toInv <- c(toInv,i)
    }# if none of these are true - marker is ok, do nothing to it!
  }
  
  ### inverting in f2: 1<->3 and 4<->5, any other case: 1<->2
  if(populationType == "f2"){ #TODO: Updated this very primitive inversion
    invertM <- newgeno[,toInv]
    invertM[which(invertM==1)] <- 3
    invertM[which(invertM==3)] <- 1
    invertM[which(invertM==5)] <- 4
    invertM[which(invertM==4)] <- 5
    newgeno[,toInv] <- invertM
  }else{ 
    newgeno[,toInv] <- 3 - newgeno[,toInv]
  }
  
  return(cbind(newgeno,oldgeno))
}

cleanGeno.internal <- function(genoCol,env,genos){
  for(envVal in unique(env)){
    incorr <- 0
    for(geno in genos){
      if(sum(which(genoCol==geno) %in% which(env==envVal)) < 4) incorr <- 1
    }
    if(incorr!=0) genoCol[which(env==envVal)] <- NA
  }
  invisible(genoCol)
}

###########################################################################################################
#                                    *** bestCorelated.internal ***
#
# DESCRIPTION:
#   subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#  which physicall map chromosome majority of markers comes
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
bestCorelated.internal <- function(cross, population, originalMap, corSDTreshold, verbose=FALSE){
  cormatrix                 <- map2mapCorrelationMatrix(cross,population,verbose)
  maximums                  <- apply(abs(cormatrix), 2, max)
  means                     <- apply(abs(cormatrix), 2, mean)
  sds                       <- apply(abs(cormatrix), 2, sd)
  selected                  <- which(maximums > (means+corSDTreshold*sds))  # Select markers that are correlated highly with more than one of the old markers
  cormatrix                 <- cormatrix[,selected]
  bestCorMarkers            <- matrix(0,length(selected),2)
  bestCorMarkers[,1]        <- apply(abs(cormatrix),2,function(r){rownames(cormatrix)[which.max(r)]})
  bestCorMarkers[,2]        <- apply(abs(cormatrix),2,function(r){rownames(cormatrix)[which.max(r[-which.max(r)])]})
  bestCorMarkers[,3]        <- maximums[selected]
  rownames(bestCorMarkers)  <- rownames(cormatrix)
  output                    <- t(apply(bestCorMarkers,1,bestCorelatedSub.internal,originalMap))
  invisible(output)
}

bestCorelatedSub.internal <- function(bestCorMarkersRow,originalMap){
  chr <- originalMap[bestCorMarkersRow[1],1]
  pos <- mean(originalMap[bestCorMarkersRow[1],2],originalMap[bestCorMarkersRow[2],2])
  invisible(c(chr,pos,bestCorMarkersRow[3]))
}

bestQTLSub.internal <- function(qtls,marker){
  cur_max <- which.max(qtls$lod[marker,])
  cur_row <- c(qtls$chr[marker,cur_max], qtls$pos[marker,cur_max], max(qtls$lod[marker,]))
  invisible(cur_row)
}

processInteractions <- function(marker,interactionType,interactionVal,threshold,flagged){
  output <- "keep"
  if(interactionVal > threshold){
    if(flagged=="warn"){
      cat("Marker:",marker,"shows significant",interactionType,"association\n")
      output <- "report"
    }else if(flagged=="remove"){
      output <- NULL
    }
  }
  invisible(output)
}

###########################################################################################################
#                                    *** bestQTL.internal ***
#
# DESCRIPTION:
#   subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#  which physicall map chromosome majority of markers comes
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
bestQTL.internal <- function(cross, population, threshold, flagged, env, verbose=FALSE, debugMode=0){
  if(is.null(population$offspring$genotypes$qtl))              stop("No qtl data in population$offspring$genotypes$qtl, run scan.qtls function first.")
  if(is.null(population$offspring$genotypes$qtl$interactions)) stop("Old version of the QTL scan detected. Re-run scan.qtls!")
  
  originalGeno          <- population$offspring$genotypes$real
  markerNames           <- markernames(cross)
  newGeno               <- pull.geno(cross)
  output                <- NULL
  
  peaksMatrix           <- getpeaks.internal(abs(population$offspring$genotypes$qtl$lod),threshold)
  rownames(peaksMatrix) <- markerNames
  
  envInt                <- 0 # nr of markers affectted by environmental interaction
  epiInt                <- 0 # nr of markers affectted by epistatic interaction
  noQTL                 <- 0 # nr of markers with no significant QTLsa
  multiQTL              <- 0 # nr of markers with multiple significant QTLs

  envMarkers <- NULL
  epiMarkers <- NULL
  
  for(marker in markerNames){
    ### is there a single significant peak in the data?
    nPeaks          <- sum(peaksMatrix[marker,]==2)
    QTLlod          <- max(population$offspring$genotypes$qtl$lod[marker,])
    envInteractions <- population$offspring$genotypes$qtl$interactions[marker,c(1,2)]
    epiInteractions <- population$offspring$genotypes$qtl$interactions[marker,3]
    affectedByEnv   <- processInteractions(marker,"environmental",max(envInteractions),(threshold/2),flagged)
    if(affectedByEnv=="report" || is.null(affectedByEnv)){
      envMarkers <- c(envMarkers,marker)
      envInt     <- envInt + 1
    }
    affectedByEpi   <- processInteractions(marker,"epistatic",max(epiInteractions),(threshold/2),flagged)
    if(affectedByEpi=="report" || is.null(affectedByEpi)){
      epiMarkers <- c(epiMarkers,marker)
      epiInt     <- epiInt + 1
    }

    curOutput <- c(NA,NA,NA)
    if(!is.null(affectedByEnv)){
      if(nPeaks==1){
        if(!is.null(affectedByEpi)){
          curOutput <- bestQTLSub.internal(population$offspring$genotypes$qtl,marker)
        }
      }else if(nPeaks < 1){
        noQTL      <- noQTL + 1
      }else{
        multiQTL   <- multiQTL + 1
      }
    }
    output <- rbind(output,curOutput)
  }
  #to have same format of the output as in bestcorrelated
  if(verbose){
    cat("\n=== Selection statistics ===\n")
    cat("Selecting from:\n\t",sum(nmar(cross)),"candidate markers.\n")
    if(flagged=="remove"){
      cat("Removed:\n")
      cat("\t",envInt,"markers showing significant association with environment.\n")
      cat("\t",epiInt,"markers influenced by an epistatic interaction.\n")
    }else if(flagged=="warn"){
      cat("\t",envInt,"markers show significant association with environment.\n")
      cat("\t",epiInt,"markers are influenced by an epistatic interaction.\n")
      cat("\nRemoved:\n")
    }
    cat("\t",noQTL,"markers showing no significant QTL.\n")
    cat("\t",multiQTL,"markers showing multiple significant QTL.\n")
  }
  rownames(output) <- markerNames
  if(any(is.na(output[,1]))) output <- output[-which(is.na(output[,1])),]
  if(any(is.na(output[,2]))) output <- output[-which(is.na(output[,2])),]
  invisible(list(output,markerNames,envMarkers,epiMarkers))
}

###########################################################################################################
#                                    *** QTLscan.internal ***
#
# DESCRIPTION:
# subfunction by Danny Arends to map QLTs modfied to work on a single phenotype
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
getpeaks.internal <- function(qtlprofiles, cutoff = 4.0){
  if(!any(qtlprofiles==Inf)) qtlprofiles[which(qtlprofiles==Inf)] <- 1000
  mmatrix <- NULL
  for(x in 1:nrow(qtlprofiles)){
    peak <- FALSE
    curmax <- 0
    curmaxindex <- 1
    marker <- 1
    maximums <- NULL
    mrow <- rep(0,ncol(qtlprofiles))
    for(ab in (qtlprofiles[x,]>cutoff | qtlprofiles[x,]<(-cutoff))){
      if(ab){
        peak <- TRUE
        if(qtlprofiles[x,marker]/abs(qtlprofiles[x,marker]) > 0){
          if(qtlprofiles[x,marker] > curmax){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- marker
          }
        }else{
          if(qtlprofiles[x,marker] < (-curmax)){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- -marker
          }
        }
        if(ncol(qtlprofiles)==marker){
          if(curmax!=0) maximums <- c(maximums,curmaxindex)
        }
      }else{
        if(curmax!=0) maximums <- c(maximums,curmaxindex)
        peak <- FALSE
        curmax <- 0
      }
      marker <- marker+1
    }
    mrow[which(qtlprofiles[x,] > cutoff)] <- 1
    mrow[which(qtlprofiles[x,] < -cutoff)] <- -1
    for(a in which(maximums>0)){ mrow[maximums[a]] <- 2 }
    for(b in which(maximums<0)){ mrow[(-maximums[b])] <- -2 }
    mmatrix <- rbind(mmatrix,mrow)
  }
  mmatrix
}

###########################################################################################################
#                                    *** map2mapCorrelationMatrix.internal ***
#
# DESCRIPTION:
#   calculating correlation matrix between genotypes inside cross object and ones from population
# OUTPUT:
#   matrix of correlations
############################################################################################################
map2mapCorrelationMatrix<- function(cross,population,verbose=FALSE){
  if(missing(cross)) stop("Please provide a cross object\n")
  if(missing(population)) stop("Please provide original genotypes\n")

  genotypes <- pull.geno(cross)
  if(verbose) cat("Calculating correlation matrix\n")
  if(!is.null(population$offspring$genotypes$real)){
    genotypesCorelationMatrix <- apply(genotypes,2,function(cgc){cor(cgc,t(population$offspring$genotypes$real),use="pair")})
    colnames(genotypesCorelationMatrix) <- colnames(genotypes)
    rownames(genotypesCorelationMatrix) <- rownames(population$offspring$genotypes$real)    
    invisible(genotypesCorelationMatrix)
  }else{
    stop("Load known genotypes into the population using add.to.population(p,genotypes,\"offspring$genotypes\")")
  }
}

###########################################################################################################
#                                           *** map2mapImage ***
#
# DESCRIPTION:
#   subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#  which physicall map chromosome majority of markers comes# 
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
map2mapImage <- function(genotypesCorelationMatrix,population,cross,corThreshold=0.5,verbose=FALSE){
  if(missing(genotypesCorelationMatrix)){
    cat("Correlation matrix not provided, calulating one")
    genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross,population,verbose)
    if(missing(cross)) stop("No object of class cross, please run either cross.denovo or enrichExistingMap\n")
    if(missing(population)) stop("Please provide a population object\n")
    check.population(population)
  }
  heatmap(genotypesCorelationMatrix,breaks = c(-1,-corThreshold,corThreshold,1),col=c("blue","white","red"),Rowv=NA)
}

