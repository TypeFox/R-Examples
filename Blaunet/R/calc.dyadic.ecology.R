calc.dyadic.ecology <-
function(blauObj, m.dist){ #splitter function
  if (m.dist == TRUE){
    blauObj$dyadic <- as.data.frame(matrix(0, nrow = 0, ncol= 8))
  }
  else{
    blauObj$dyadic <- as.data.frame(matrix(0, nrow = 0, ncol= 7))
  }

  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]

  for(ecologyId in uniqueEcologies) {
    ecologyRows <- which(blauObj$ids[,2] == ecologyId)
    
    miniBlau <- splittify(blauObj, ecologyId, ecologyRows)

    miniBlau <- calc.dyadic(miniBlau, m.dist)

    blauObj$dyadic <- rbind(blauObj$dyadic, miniBlau$dyadic)

  }
  return(blauObj)
}
