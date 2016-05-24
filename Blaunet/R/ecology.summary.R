ecology.summary <-
function(blauObj, percent = FALSE){
  uniqueEcologies <- unique(blauObj$ids[,2])
  uniqueEcologies <- uniqueEcologies[!is.na(uniqueEcologies)]
  sum.ecology <- matrix(0, nrow = 0, ncol = (ncol(blauObj$memberships) + 2)) #ecology name, niche name, num in org, num in niche, num in org but not niche, num not in any other niche, overlap 2+ other niches ##no boundaries, too cluttered
  colnames(sum.ecology) <- c("Ecology", "Org/Niche", colnames(blauObj$memberships))
  if (length(uniqueEcologies)==1) {
    focalNiches <- blauObj$isInNiche
    sum.mat <- t(as.matrix(focalNiches)) %*% as.matrix(focalNiches)
    orig.diag <- diag(sum.mat)
    #how many individuals are exclusively in a niche?
    #get this by summing up all isInNiche rows with sum = 1
    mat.diagonal <- as.matrix(rep(0, ncol(focalNiches)))
    for (node in 1:nrow(focalNiches)){
      if(sum(focalNiches[node,]) == 1){
        mat.diagonal <- mat.diagonal + focalNiches[node,]
      }
    }
    #manually replace because diag function gets confused
    for (elem in 1:nrow(mat.diagonal)){
      sum.mat[elem,elem] <- mat.diagonal[elem,1]
    }
    if (percent == TRUE){
      sum.mat <- sum.mat / orig.diag
    }
    sum.ecology <- rbind(sum.ecology, cbind(cbind(rep(blauObj$ids[1,2], ncol(blauObj$memberships)), colnames(blauObj$memberships)), sum.mat))
  } else {
    for (ecologyId in uniqueEcologies){
      ecologyRows <- which(blauObj$ids[,2] == ecologyId)
      focalNiches <- blauObj$isInNiche[ecologyRows, 1:(ncol(blauObj$isInNiche)-1)]
      sum.mat <- t(as.matrix(focalNiches)) %*% as.matrix(focalNiches)
      orig.diag <- diag(sum.mat)
      #how many individuals are exclusively in a niche?
      #get this by summing up all isInNiche rows with sum = 1
      mat.diagonal <- as.matrix(rep(0, ncol(focalNiches)))
      for (node in 1:nrow(focalNiches)){
        if(sum(focalNiches[node,]) == 1){
          mat.diagonal <- mat.diagonal + focalNiches[node,]
        }
      }
      #manually replace because diag function gets confused
      for (elem in 1:ncol(mat.diagonal)){
        sum.mat[elem,elem] <- mat.diagonal[1,elem]
      }
      if (percent == TRUE){
        sum.mat <- sum.mat / orig.diag
      }
      sum.ecology <- rbind(sum.ecology, cbind(cbind(rep(ecologyId, ncol(blauObj$memberships)), colnames(blauObj$memberships)), sum.mat))
    }
  }
  rownames(sum.ecology) <- NULL
  return(sum.ecology)
}
