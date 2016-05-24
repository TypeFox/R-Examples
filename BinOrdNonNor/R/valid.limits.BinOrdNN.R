valid.limits.BinOrdNN <-
function(plist, skew.vec, kurto.vec, no.bin, no.ord, no.NN) 
{
  no.binord <- no.bin + no.ord
  if (no.binord>0) validate.plist(plist, no.binord)
  minmat <- maxmat <- diag(length(plist) + no.NN)
  
  for (r in 2:nrow(minmat)) {
    for (c in 1:(r - 1)) {
      if (r <= length(plist)) {
        minmax <- Limit_forOO(plist[[r]], plist[[c]])
      }else if (r > length(plist) & c > length(plist)) {
        minmax <- Limit_forNN(skew.vec[c(r-no.binord,c-no.binord)],
                              kurto.vec[c(r-no.binord,c-no.binord)])
      }else if (r > length(plist) & c <= length(plist)) {
        minmax = Limit_forONN(plist[[c]], skew.vec[r-no.binord], kurto.vec[r-no.binord])
      }
      minmat[r, c] <- minmax[1]
      maxmat[r, c] <- minmax[2]
      rm(minmax)
    }
  }
  minmat <- minmat + t(minmat)
  diag(minmat) <- 1
  maxmat <- maxmat + t(maxmat)
  diag(maxmat) <- 1
  return(list(lower = minmat, upper = maxmat))
}
