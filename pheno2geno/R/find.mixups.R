#
# find.mixups.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Jan, 2012
# Contains: cross.saturate, rearrangeMarkers, bestCorelated.internal, map2mapImage
#           map2mapCorrelationMatrix.internal 
#

# find.mixups
#
# DESCRIPTION:
#  Saturate existing genetic map adding markers derived from gene expression
# OUTPUT:
#  An object of class cross
#
find.mixups <- function(population,map=c("genetic","physical"),n.qtls=50,threshold=15,verbose=FALSE){
  s <- proc.time()
  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using add.to.population function\n")
  }
  if(!is.finite(threshold)){
    stop("Threshold must be a finite, numeric value.\n")
  }
  if(!is.finite(n.qtls)){
    stop("n.qtls must be a finite, numeric value.\n")
  }
  if(n.qtls<0 || n.qtls>nrow(population$offspring$phenotypes)){
    stop("Value of n.qtls is too high or too low.\n")
  }
  map <- match.arg(map)
  if(map=="genetic"){
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$genetic))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$genetic <- population$maps$genetic[rownames(population$offspring$genotypes$real),]
      if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
    }
    ### THAT's just ugly trick to make saving the cross with a lot fo phenotypes faster.
    population10pheno <- population
    population10pheno$offspring$phenotypes <- population10pheno$offspring$phenotypes[1:10,]
    ### creation of the cross
    tryCatch({
      aa <- tempfile()
      sink(aa)
      returncross <- genotypesToCross.internal(population10pheno,"real","map_genetic")
      returncross$pheno <- t(population$offspring$phenotypes)
    },
    error= function(err){
      stop(paste("ERROR in find.mixups while creating cross:  ",err))
      sink()            # sink if errored -> otherwise everything is sinked into aa file
      # file is not removed -> contains output that may help with debugging
    },
    finally={
      sink()
      file.remove(aa) # no error -> close sink and remove unneeded file
    })
  }else{
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$physical))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$physical <- population$maps$physical[rownames(population$offspring$genotypes$real),]
      if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
    }
    #for faster creation of cross
    population10pheno <- population
    population10pheno$offspring$phenotypes <- population10pheno$offspring$phenotypes[1:10,]
    ### creation of the cross
    tryCatch({
      aa <- tempfile()
      sink(aa)
      returncross <- genotypesToCross.internal(population10pheno,"real","map_physical")
      returncross$pheno <- t(population$offspring$phenotypes)
    },
    error= function(err){
      stop(paste("ERROR in find.mixups while creating cross:  ",err))
      sink()            # sink if errored -> otherwise everything is sinked into aa file
      # file is not removed -> contains output that may help with debugging
    },
    finally={
      sink()
      file.remove(aa) # no error -> close sink and remove unneeded file
    })
  }
  returncross <- calc.genoprob(returncross)
  qtls_found <- 0
  qtls <- NULL
  markers <- rownames(population$offspring$phenotypes)
  scores <- vector(mode="numeric",length=ncol(population$offspring$phenotypes))
  names(scores) <- colnames(population$offspring$phenotypes)
  done <- 0
  if(verbose) cat("--- starting the QTL analysis ---\n")
  while(qtls_found<n.qtls){
    phenotype <- round(runif(1,1,nrow(population$offspring$phenotypes)))
    cur_phenotype <- matrix(scanone(returncross,pheno.col=phenotype,method="hk")[,3],1,nrow(population$offspring$genotypes$real))
    cur_peaks <- getpeaks.internal(abs(cur_phenotype),threshold)
    done <- done+1
    if(any(cur_peaks==2)){
      peakLocations <- which(cur_peaks==2)
      qtls_found <- qtls_found + length(peakLocations)
      if(verbose) cat(qtls_found,"qtls found, phenotype:",phenotype,"marker:",peakLocations,"\n")
      old_names <- names(qtls)
      qtls <- c(qtls,rep(phenotype,length(peakLocations)))
      names(qtls) <- c(old_names,peakLocations)
    }
  }
  newGeno <- population$offspring$genotypes$real
  for(j in 1:qtls_found){
    group_a <- which(population$offspring$genotypes$real[as.numeric(names(qtls))[j],]==1)
    group_b <- which(population$offspring$genotypes$real[as.numeric(names(qtls))[j],]==2)
    newGeno_ <- scoreMixups.internal(group_a,group_b,scores,qtls_found,population$offspring$phenotypes[qtls[j],],newGeno,as.numeric(names(qtls))[j])
    newGeno <- newGeno_[[1]]
    scores <- newGeno_[[2]]
  }
  if(verbose){
    cat("Numbers of the phenotypes scanned:",done,"\n")
    if(any(scores>50)){
      flagged <- which(scores>50)
      cat("Found",length(flagged),"possible mix-ups:\n")
      for(flag in flagged){      
        cat(names(scores)[flag],":",scores[flag],"% flagged\n")
      }
    }
  }
  e <- proc.time()
  if(verbose) cat("Function runtime:",(e-s)[3],"s \n")
  invisible(list(newGeno,scores))
}

scoreMixups.internal <- function(group_a,group_b,scores,qtls_found,curRow,genotypes,n.curMarker){
  newRow <- genotypes[n.curMarker,]
  meanGroupA <- mean(curRow[group_a])
  meanGroupB <- mean(curRow[group_b])
  rowMean <- mean(curRow)
  increase <- 1/qtls_found*100
  if(meanGroupA>meanGroupB){
    if(any(curRow[group_a]<rowMean)){
      positions <- names(curRow[group_a])[which(curRow[group_a]<rowMean)]
      scores[positions] <- scores[positions]+increase
      newRow[positions] <- 2
    }
    if(any(curRow[group_b]>rowMean)){
      positions <- names(curRow[group_b])[which(curRow[group_b]>rowMean)]
      scores[positions] <- scores[positions]+increase
      newRow[positions] <- 1
    }
  }else{
    if(any(curRow[group_a]>rowMean)){
      positions <- names(curRow[group_a])[which(curRow[group_a]>rowMean)]
      scores[positions] <- scores[positions]+increase
      newRow[positions] <- 2
    }
    if(any(curRow[group_b]<rowMean)){
      positions <- names(curRow[group_b])[which(curRow[group_b]<rowMean)]
      scores[positions] <- scores[positions]+increase
      newRow[positions] <- 1
    }
  }
  genotypes[n.curMarker,] <- newRow
  invisible(list(genotypes,scores))
}
