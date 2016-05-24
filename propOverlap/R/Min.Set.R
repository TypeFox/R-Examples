Min.Set <- function(GMask, Scores){
  Grand.Mask <- apply(GMask, 2, function(x){any(x==1)})
  Grand.Min  <- rep(FALSE, length(Grand.Mask))
  Min.Genes  <- c()
  k = 0
  while (any(Grand.Min != Grand.Mask)){
    k = k + 1 
    Sum.ones     <- GMask %*% rep(1, ncol(GMask))
    Temp.Set     <- which(Sum.ones == max(Sum.ones))
    Min.Genes[k] <- Temp.Set[which.min(Scores[Temp.Set])]
    Grand.Min    <- Grand.Min | (GMask[Min.Genes[k],]==1)
    GMask[, Grand.Min] = 0
  }
  return(list(Min.Subset = Min.Genes, Covered.Obs = which(Grand.Min)))
}