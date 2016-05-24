calc.niches.ecology <-
function(blauObj, uniqueEcologies, dev.range){
  
  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]
  
  blauObj$isInNiche <- matrix(0, nrow = nrow(blauObj$dimensions), ncol = (ncol(blauObj$memberships) + 1)) #extra column for ecology names
  colnames(blauObj$isInNiche) <- c(vapply(colnames(blauObj$memberships), function(x) paste(x, "niche", sep="_"), "a"), 'ecologyNames')
  rownames(blauObj$isInNiche) <- blauObj$ids[,1]
  
  for(ecologyId in uniqueEcologies){ #iterate through each ecology: all of the calculations for the ecology happen here and they are appended to $isInNiche, $topbounds, and $lowbounds
    ecologyRows <- which(blauObj$ids[,2] == ecologyId) #pull out ROW identifiers for each row in the ecology
    
    miniBlau <- splittify(blauObj, ecologyId, ecologyRows)

    miniBlau <- calc.niches(miniBlau, dev.range) #memberships, dimensions, primaryMemberships are used by niches
    
    blauObj$isInNiche[ecologyRows,] <- cbind(miniBlau$isInNiche, (rep(ecologyId, nrow(miniBlau$isInNiche))))
    
    topbounds <- cbind(as.data.frame(miniBlau$topbounds), as.data.frame(rep(ecologyId, nrow(miniBlau$topbounds)))) #temp object to add ecology names 
    lowbounds <- cbind(as.data.frame(miniBlau$lowbounds), as.data.frame(rep(ecologyId, nrow(miniBlau$lowbounds)))) #temp object to add ecology names
    
    colnames(topbounds) <- c(colnames(blauObj$dimensions), 'ecologyNames')
    rownames(topbounds) <- colnames(blauObj$memberships)
    colnames(lowbounds) <- c(colnames(blauObj$dimensions), 'ecologyNames')
    rownames(lowbounds) <- colnames(blauObj$memberships)
    
    blauObj$topbounds <- rbind(blauObj$topbounds, topbounds) #add it to the bottom
    blauObj$lowbounds <- rbind(blauObj$lowbounds, lowbounds) #add it to the bottom
    
  }

  tempData <- blauObj$isInNiche[,which(colnames(blauObj$isInNiche) != 'ecologyNames')]
  class(tempData) <- 'numeric'

  blauObj$isInNiche <- cbind(as.data.frame(tempData), blauObj$isInNiche[,which(colnames(blauObj$isInNiche) == 'ecologyNames')])

  colnames(blauObj$isInNiche)[ncol(blauObj$isInNiche)] <- 'ecologyNames'

  return(blauObj)

}
