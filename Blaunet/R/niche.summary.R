niche.summary <-
function(blauObj){

  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]

  sum.niche <- matrix(0, nrow = length(uniqueEcologies)*ncol(blauObj$memberships), ncol = 7) #ecology name, niche name, num in org, num in niche, num in org but not niche, num not in any other niche, overlap 2+ other niches ##no boundaries, too cluttered

  colnames(sum.niche) <- c("Ecology", "Org/Niche", "OrgMem", "NicheMem", "NicheExc", "NicheOvr", "MemExc")

  rowCount <- 1

  if (length(uniqueEcologies)==1) {
    nicheNum <- 1
    focalMemberships <- blauObj$memberships
    focalNiches <- blauObj$isInNiche
    for (colCyc in 1:ncol(focalMemberships)){
      #in org but not in niche
      #basically find the 1's
      diff <- focalMemberships[, nicheNum] - focalNiches[, nicheNum]
      numOutside <- length(which(diff == 1))
      #in focal niche but in no other niche
      nicheExcl <- sum(apply(focalNiches, 1, function(x) ifelse(sum(x) == 1 && x[nicheNum] == 1, 1, 0)), na.rm = TRUE)
      #overlaps with at least 2 niches
      numNonExclusive <- sum(focalNiches[, nicheNum], na.rm = TRUE) - nicheExcl
      sum.niche[rowCount,] <- c(blauObj$ids[1,2], colnames(blauObj$memberships)[nicheNum], sum(focalMemberships[, nicheNum], na.rm = TRUE), sum(focalNiches[, nicheNum], na.rm = TRUE), nicheExcl, numNonExclusive, numOutside)
      rowCount <- rowCount + 1
      nicheNum <- nicheNum + 1
    }
  } else {
    for (ecologyId in uniqueEcologies){
      nicheNum <- 1
      ecologyRows <- which(blauObj$ids[,2] == ecologyId)
      focalMemberships <- blauObj$memberships[ecologyRows, ]
      focalNiches <- blauObj$isInNiche[ecologyRows, 1:(ncol(blauObj$isInNiche)-1)] #exclude last column, which is ecology index
      for (colCyc in 1:ncol(focalMemberships)){
        #in org but not in niche
        #basically find the 1's
        diff <- focalMemberships[, nicheNum] - focalNiches[, nicheNum]
        numOutside <- length(which(diff == 1))
        #in focal niche but in no other niche
        nicheExcl <- sum(apply(focalNiches, 1, function(x) ifelse(sum(x) == 1 && x[nicheNum] == 1, 1, 0)), na.rm = TRUE)
        #overlaps with at least 2 niches
        numNonExclusive <- sum(focalNiches[, nicheNum], na.rm = TRUE) - nicheExcl
        sum.niche[rowCount,] <- c(ecologyId, colnames(blauObj$memberships)[nicheNum], sum(focalMemberships[, nicheNum], na.rm = TRUE), sum(focalNiches[, nicheNum], na.rm = TRUE), nicheExcl, numNonExclusive, numOutside)
        rowCount <- rowCount + 1
        nicheNum <- nicheNum + 1
      }
    }
  }
  return(sum.niche)
}
