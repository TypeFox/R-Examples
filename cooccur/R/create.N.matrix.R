create.N.matrix <-
function(mat){
  nspp <- nrow(mat)
  Nmat <- matrix(nrow = nspp,ncol = nspp)
  row <- 0
  for (spp in 1:nspp){
    if (spp < nspp){
      for (spp_next in (spp + 1):nspp){
        Nmat[spp,spp_next] <- sum(mat[spp,]*mat[spp_next,])
        Nmat[spp_next,spp] <- sum(mat[spp,]*mat[spp_next,])
      }  
    }
  }  
  return(Nmat)
}
