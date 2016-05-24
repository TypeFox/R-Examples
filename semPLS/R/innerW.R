# Constructs the adjacency matrix D of the structural model.
innerW <-
function(strucmod, latent){
  D <- matrix(0, ncol=length(latent), nrow=length(latent))
  colnames(D) <- latent             # names of LVs
  rownames(D) <- latent             # names of LVs
  for(i in 1:nrow(strucmod)){
    D[which(latent==strucmod[i,1]), which(latent==strucmod[i,2])] <- 1
  }
  return(D)
}
