predict.bayescomm <- function(object, newdata, ...) {
  # Return an array of sample probabilities at new sites.  
  # Rows are test sites; columns are species; slices are coefficient samples
 
  bindSpeciesCoefficients <- function(object) {
    # Helper function to bind the coefficient lists in B into an array.
    
    pre.binding <- lapply(object$trace$B, function(x) {
      dim(x) <- c(dim(x), 1)
      x
    })
    
    out <- do.call(abind, pre.binding)
    colnames(out) <- colnames(object$trace$B[[1]])
    
    out
  } 
  
  
  # Probably only works with the "full" model type and no `covlist` or
  # `condition` specified.
  
  # Haven't played with any cases where mu isn't null
  if(!is.null(object$other$mu)){
    stop("predictions are not supported for non-null mu")
  }
  
  X <- cbind(intercept = 1, newdata)
  B <- bindSpeciesCoefficients(object)
  R <- object$trace$R
  n.species <- dim(B)[3]
  
  predictions <- array(
    NA, 
    dim = c(nrow(X), dim(B)[3], nrow(B)), 
    dimnames = list(row.names(X), dimnames(B)[[3]], NULL)
  )
  
  # Fill in predictions slice by slice
  for (i in 1:nrow(B)) {
    Sigma <- matrix(0, nrow = n.species, ncol = n.species)
    Sigma[upper.tri(Sigma)] <- R[i, ]  # Fill in upper triangle?
    Sigma <- Sigma + t(Sigma)  # Fill in lower triangle?
    diag(Sigma) <- 1  # Diagonal equals 1 in multivariate probit model
    
    Z <- rmvnorm(n = nrow(X), mean = rep(0, n.species), sigma = Sigma)
    predictions[, , i] <- pnorm(X %*% B[i, , ] + Z)
  }
  
  predictions
}