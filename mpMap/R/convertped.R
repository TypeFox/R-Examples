convertped <- function(pedigree)
{
  # convert a pedigree to numeric
  ped <- apply(pedigree, 2, as.character)
  mid <- as.numeric(match(ped[,2], ped[,1]))
  fid <- as.numeric(match(ped[,3], ped[,1]))

  # Any parents not in the pedigree are unknown
  mid[is.na(mid)] <- 0
  fid[is.na(fid)] <- 0

  pednum <- cbind(1:nrow(ped), mid, fid, as.numeric(ped[,4]))

#  pednum <- as.data.frame(pednum)
  colnames(pednum) <- colnames(pedigree)

  return(pednum)
}

