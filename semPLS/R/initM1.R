# Used in 'read.splsm':
# Constructs the initial outer weights.
initM1 <-
function(model){
  measuremod <- model$measuremod
  latent <- model$latent                              # names of LVs
  mf <- model$manifest                                # names of MVs
  M <- matrix(0, ncol=length(latent), nrow=length(mf))
  colnames(M) <- latent
  rownames(M) <- mf
  for(i in 1:nrow(measuremod)){
    if (measuremod[i,1] %in% mf)
      M[which(mf==measuremod[i,1]), which(latent==measuremod[i,2])] <- 1
    if (measuremod[i,2] %in% mf)
      M[which(mf==measuremod[i,2]), which(latent==measuremod[i,1])] <- 1    
  }
  return(M)
}
