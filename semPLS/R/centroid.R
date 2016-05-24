#' A function for the centroid fitting scheme for the inner weights within sempls. The centroid scheme only uses the sign of the correlation between latent variables connected with each other.
centroid <-
function(model, fscores, pairwise, method){
  method <- "pearson" ## for other methods: convergence problems!
  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  D <- model$D
  C <- (D + t(D))
  innerW <- C
  innerW[C == 1] <- sign(cor(fscores, use=use, method=method))[C == 1]
  return(innerW)
}
