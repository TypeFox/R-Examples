# Initial estimation of factor scores
step1 <-
function(model, data, sum1, pairwise, method, ...){
  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  M <- model$M
  if(sum1) M <- apply(M, 2, sum1)
  if(pairwise){
    blocks <- model$blocks
    nl <- length(model$latent) # number of LVs
    #root <- vector(mode="list", length=nl)
    #names(root) <- model$latent
    Latent <- matrix(NA, nrow=nrow(data), ncol=nl)
    colnames(Latent) <- model$latent
    for(i in model$latent){
      if(length(blocks[[i]])==1){
          #root[[i]] <- 1
          Latent[,i] <- as.matrix(data[ , blocks[[i]] ])
          next
      }
      mf <- as.matrix(data[ , blocks[[i]]])        # MVs in i-th LVs block
      #root[[i]] <- solve(chol(cor(mf,y=NULL, use, method)))
      w <- as.matrix(M[blocks[[i]], i])
      #w <- root[[i]] %*% w/norm(w, "F")
      M[blocks[[i]], i] <- w
      Latent[,i] <- mf %*% w
    }
    Latent <- scale(Latent)
  }
  else {Latent <- scale(as.matrix(data) %*% M)}
  # the attributes for the scale are meaningless
  #attributes(Latent)[c(3,4)] <- NULL
  return(list(latent=Latent, outerW=M))
}
