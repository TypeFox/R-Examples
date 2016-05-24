dyadic <-
function(blauObj, dev.range = 1.5, ecologies.off=FALSE, m.dist = TRUE) {

  if (ecologies.off == TRUE){
    blauObj <- niches(blauObj, dev.range, ecologies.off)
  }
  
  uniqueEcologies <- unique(blauObj$ids[,2])
  if(length(uniqueEcologies) == 1 || ecologies.off == TRUE){ #if there's only one ecology and ,we don't have isInNiche

    if(is.null(blauObj$isInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }

    blauObj <- calc.dyadic(blauObj, m.dist)

  }

  else if(length(uniqueEcologies) > 1){
    if(is.null(blauObj$IsInNiche)){
      blauObj <- niches(blauObj, dev.range)
    }
    blauObj <- calc.dyadic.ecology(blauObj, m.dist)
  }

  if (m.dist == TRUE){
    colnames(blauObj$dyadic) <- c('Ego', 'Alter', 'CoNicher', 'CoOutsider', 'Straddler', 'Spanner', 'EucDist', 'MahalanobisDist')
  }
  else{
    colnames(blauObj$dyadic) <- c('Ego', 'Alter', 'CoNicher', 'CoOutsider', 'Straddler', 'Spanner', 'EucDist')
  }

  return(blauObj)
}
