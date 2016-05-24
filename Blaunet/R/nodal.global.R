nodal.global <-
function(blauObj, dev.range=1.5, ecologies.off = FALSE){
  
  if (ecologies.off == TRUE){
    blauObj <- niches(blauObj, dev.range, ecologies.off)
  }

  uniqueEcologies <- unique(blauObj$ids[,2])
  if(length(uniqueEcologies) == 1 || ecologies.off == TRUE){ #if there's only one ecology and ,we don't have isInNiche

    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }

    blauObj <- calc.nodal(blauObj, mode = "global") #has isInNiches now; does a bunch of stuff if a primaryMembership is specified
  }

  else if(length(uniqueEcologies) > 1) { #if there's more than one ecology, we need to split up primary membership, weights, dimensions and memberships
    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }    
    blauObj <- calc.nodal.ecology(blauObj, uniqueEcologies, mode = "global")
  }
  colnames(blauObj$nodalGlobal) <- c("TotalOrgs", "Nicher", "NicheList")
  
  return(blauObj)
}
