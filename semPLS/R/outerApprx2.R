# Outer Approximation (step 3 within PLS-Algorithm)
# replaces: outerApprx() [better implementation]
# Calculates the new outer weights.
# uses: sum1 - function to normalize the weights to sum up to 1.
# last modified: 05.05.2011 (Armin Monecke)
outerApprx2 <-
function(Latent, data, model, sum1, pairwise, method){
  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  blocks <- model$blocks
  N <- nrow(data)
  nl <- ncol(Latent)                      # number of latent variables
  W <- model$M
  for (i in model$latent){
    if(length(blocks[[i]])==1) next
    mf <- as.matrix(subset(data, select=blocks[[i]]))
    fscores <- as.matrix(Latent[,i])

    # only for Mode "B": transform the MVs of a block
    if (attr(blocks[[i]], "mode") == "B") {
      S <- cor(mf, mf, use, method)
      T <- solve(chol(S))
      mf <-  t(t(T) %*% t(mf)) 
    }
    
    # the same for mode "A" and "B"
    W[blocks[[i]],i] <- cor(fscores, mf, use, method)

    # only for Mode "B": retransform the weights according to the MVs
    if (attr(blocks[[i]], "mode") == "B") {
      W[blocks[[i]],i] <- T %*% W[blocks[[i]],i]
    }
  }

  ## Normalize weights to colwise sum up to 1?
  if(sum1==TRUE){
     W <- apply(W, 2, sum1)
  }
  return(W)
}
