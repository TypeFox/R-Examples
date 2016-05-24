# factorial weighting scheme
factorial <- function(model, fscores, pairwise, method){
  method <- "pearson" ## for other methods: convergence problems!
  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  D <- model$D
  C <- D + t(D)
  innerW <- C
  innerW[C == 1] <- cor(fscores, use=use, method=method)[C == 1]
  return(innerW)
}
