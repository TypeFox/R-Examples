# path weighting scheme
# major changes from revision 16 to 17
pathWeighting <-
function(model, fscores, pairwise, method){
  method <- "pearson" ## for other methods: convergence problems!
  ifelse(pairwise, use <- "pairwise.complete.obs",
                   use <- "everything")
  D <- model$D
  latent <- model$latent
  E <- D - t(D)
  pred <- predecessors(model)
  # calculating the inner weights
  innerW <- E
  for (i in latent){
    if(length(pred[[i]])==0) next
    else if (length(pred[[i]])==1){
        innerW[pred[[i]], i] <- cor(fscores[,pred[[i]]], fscores[,i],
                                    use=use, method=method)
    }
    innerW[pred[[i]], i] <- solve(cor(as.matrix(fscores[,pred[[i]]])
                                      , use=use, method=method)) %*%
                                  cor(fscores[,pred[[i]]], fscores[,i], use=use, method=method)

  }

  innerW[E == 0] <- 0
  innerW[E == -1] <- cor(as.matrix(fscores[, latent]), use=use, method=method)[E == -1]
  # return the matrix of inner weights
  return(innerW)
}
