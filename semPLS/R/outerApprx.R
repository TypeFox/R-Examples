# Outer Approximation (step 3 within PLS-Algorithm)
# Calculates the new outer weights.
# uses: sum1 - function to normalize the weights to sum up to 1.
# last modified: 21.10.2010 (Armin Monecke)
outerApprx <-
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
    ## Mode A: reflective
    if (attr(blocks[[i]], "mode") == "A") {
      #w <- cor(fscores, mf, use, method)                 # new
      W[blocks[[i]],i] <- cor(fscores, mf, use, method) # old
    }
    ## Mode B: formative
    else if (attr(blocks[[i]], "mode") == "B") {
      #w <- t(solve(cor(mf, mf, use, method)) %*%                 # new
      #     cor(mf, fscores, use, method))                       #new
      W[blocks[[i]],i] <- solve(cor(mf, mf, use, method)) %*% # old
                          cor(mf, fscores, use, method)       # old
    }
    # changed: 21.10.2010; Idea: constrain w_i'S_{ii}w_i=1
    # see T.K. Dijkstra, Latent Variables and Indices: Herman Wold's
    #     Basic Design and Partial Least Squares, Handbook of Partial
    #     Least Squares, p.32-33.
    # Note: cholesky decomposition for a block does not change -> put outside loop!!!
    #w <- solve(chol(cor(mf * rep(w, each=N), y=NULL, use, method))) %*%
    #     t(w/norm(w, "F"))                                                # new
    #W[blocks[[i]],i] <- w                                                # new
    ## comment: did not work!!!
  }

  ## Normalize weights to colwise sum up to 1?
  if(sum1==TRUE){
     W <- apply(W, 2, sum1)
  }
  return(W)
}
