niches <-
function(blauObj, dev.range = 1.5, ecologies.off = FALSE){
  uniqueEcologies <-  unique(blauObj$ids[,2])

  if (length(uniqueEcologies) == 1 || ecologies.off == TRUE){
    blauObj <- calc.niches(blauObj, dev.range)
    rownames(blauObj$isInNiche) <- rownames(blauObj$memberships)
  }
  else if (length(uniqueEcologies)) {
    blauObj <- calc.niches.ecology(blauObj, uniqueEcologies, dev.range)
  }

  presentCases <- which(complete.cases(blauObj$dimensions))
  
  blauObj <- getPresentCases(blauObj, presentCases)

  return(blauObj)
}
