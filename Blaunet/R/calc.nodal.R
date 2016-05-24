calc.nodal <-
function(blauObj, mode){ 

  #requires focalNiche (primMem) specification
  if (mode == 'local'){

    #initialize
    blauObj$nodalLocal <- as.data.frame(matrix(0, nrow = nrow(blauObj$memberships), ncol= 3))

    rownames(blauObj$nodalLocal) <- rownames(blauObj$isInNiche)

    #in focal niche
    blauObj$nodalLocal[,1] <- blauObj$isInNiche[, blauObj$primMemCol]

    #total number of niches individual is in
    blauObj$nodalLocal[,2] <- matrix(apply(blauObj$isInNiche, 1, function(x) sum(x, na.rm=TRUE)), ncol = 1, byrow = TRUE)

    #if individual is in primary org but outside of primary niche
    for(nodeCyc in 1:length(blauObj$isInNiche[, blauObj$primMemCol])){
      if (!is.na(blauObj$isInNiche[nodeCyc, blauObj$primMemCol]) && !is.na(blauObj$memberships[nodeCyc, blauObj$primMemCol])){
        if (blauObj$isInNiche[nodeCyc, blauObj$primMemCol] == 0 && blauObj$memberships[nodeCyc, blauObj$primMemCol] == 1){
          blauObj$nodalLocal[nodeCyc, 3] <- 1
        }
      }
    }
    return(blauObj)
  }

  #does not require focalNiche(primMem)
  else if (mode == 'global'){

    #number of organizations individual is in
    orgs <- matrix(apply(blauObj$memberships, 1, function(x) sum(x, na.rm = TRUE)), ncol = 1, byrow = TRUE)

    #number of niches individual is in
    niches <- matrix(apply(blauObj$isInNiche, 1, function(x) c(sum(x, na.rm = TRUE), c(paste(which(x == 1), collapse=' ')))), ncol = 2, byrow = TRUE)

    blauObj$nodalGlobal <- cbind(orgs, niches)

    rownames(blauObj$nodalGlobal) <- rownames(blauObj$isInNiche)

    return(blauObj)
  }
}
