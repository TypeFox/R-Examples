CR_calcLD <-
function(finals, founders, pedigree, pairs, rmat)
{
  naigen <- nai(pedigree)  

  npairs <- nrow(pairs)
  nfounders <- nrow(founders)
  nped <- nrow(pedigree)
  nfinals <- nrow(finals)
  nmrk <- ncol(founders)
  finals[,1] <- as.numeric(as.character(finals[,1]))
 
  finals[is.na(finals)] <- -1

  finvec <- as.vector(t(finals))
  fouvec <- as.vector(t(founders))

  test <- .C("calcLD", as.integer(finvec), as.integer(fouvec), as.integer(pedigree[,1]), as.integer(pedigree[,2]), as.integer(pedigree[,3]),  as.integer(pairs[,1]), as.integer(pairs[,2]), as.integer(nfinals), as.integer(nfounders), as.integer(nmrk), as.integer(npairs), as.integer(nped), as.integer(naigen), as.double(lowerTriangle(rmat, diag=TRUE)), ldW=double(length=npairs), ldLew=double(length=npairs), lddelta=double(length=npairs), ldr2=double(length=npairs), PACKAGE="mpMap")
  
  ldtable <- cbind(test$ldW, test$ldLew, test$lddelta, test$ldr2)
  colnames(ldtable) <- c("W", "LewontinD", "Delta2", "r2")

  return(ldtable)
}

