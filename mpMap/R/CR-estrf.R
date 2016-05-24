CR_estrf <-
function(finals, founders, pedigree, pairs, r)
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

  test <- .C("rfhaps", as.integer(finvec), as.integer(fouvec), as.integer(pedigree[,1]), as.integer(pedigree[,2]), as.integer(pedigree[,3]),  as.integer(pairs[,1]), as.integer(pairs[,2]), as.integer(nfinals), as.integer(nfounders), as.integer(nmrk), as.integer(npairs), as.integer(nped), as.integer(naigen), as.double(r), as.integer(length(r)), pairth=double(length=length(r)*npairs), PACKAGE="mpMap")

  return(matrix(data=test$pairth, nrow=npairs, ncol=length(r), byrow=TRUE))
}

