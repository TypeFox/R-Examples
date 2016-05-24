# Outer estimation of the factor scores
step4 <-
function(data, outerW, model, pairwise){
  blocks <- model$blocks
  if(pairwise){
    Latent <- matrix(NA, nrow=nrow(data), ncol=length(model$latent)) # factor scores
    colnames(Latent) <- model$latent
    for(i in model$latent){
      mf <- as.matrix(data[ , blocks[[i]] ])
      #Latent[,i] <- mf %*% as.matrix(outerW[blocks[[i]], i])
      Latent[,i] <- mf %*% outerW[blocks[[i]], i, drop=FALSE]
    }
    Latent <- scale(Latent)
  }
  else {Latent <- scale(as.matrix(data) %*% outerW)}  # old
  # the attributes for the scale are meaningless
  # No, they are meaningfull: w'Sw=1
  #attributes(Latent)[c(3,4)] <- NULL
  
  ## Alternatively: without scaling in each iteration
  # else {Latent <- scale(as.matrix(data) %*% outerW, center = TRUE, scale = FALSE)
  #     attr(Latent, "scaled:scale") <- 1}
  
  return(Latent)
}
