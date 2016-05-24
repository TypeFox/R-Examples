calc.nodal.ecology <-
function(blauObj, uniqueEcologies, mode){

  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]
  
  if (mode == 'local'){
    blauObj$nodalLocal <- matrix(0, nrow = 0, ncol = 2)

    for(ecologyId in uniqueEcologies){

      ecologyRows <- which(blauObj$ids[,2] == ecologyId)
      miniBlau <- splittify(blauObj, ecologyId, ecologyRows)
      
      miniBlau <- calc.nodal(miniBlau, mode) #memberships, weights, dimensions, primaryMemberships are used by niches
      blauObj$nodalLocal <- rbind(blauObj$nodalLocal, miniBlau$nodalLocal)
    }
  }

  else if (mode == 'global'){
    blauObj$nodalGlobal <- matrix(0, nrow = 0, ncol = 3)

    for(ecologyId in uniqueEcologies){

      ecologyRows <- which(blauObj$ids[,2] == ecologyId)
      miniBlau <- splittify(blauObj, ecologyId, ecologyRows)
      
      miniBlau <- calc.nodal(miniBlau, mode) #memberships, weights, dimensions, primaryMemberships are used by niches
      blauObj$nodalGlobal <- rbind(blauObj$nodalGlobal, miniBlau$nodalGlobal)
    }
  }

  if (mode == 'network'){
    blauObj$nodalNetwork <- matrix(0, nrow = 0, ncol = 2)

    for(ecologyId in uniqueEcologies){

      ecologyRows <- which(blauObj$ids[,2] == ecologyId)
      miniBlau <- splittify(blauObj, ecologyId, ecologyRows)
      
      miniBlau <- calc.nodal.network(miniBlau) #memberships, weights, dimensions, primaryMemberships are used by niches
      blauObj$nodalNetwork <- rbind(blauObj$nodalNetwork, miniBlau$nodalNetwork)
    }
  }
  return(blauObj)
}
