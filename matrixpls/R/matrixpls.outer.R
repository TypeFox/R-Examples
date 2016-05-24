# =========== Outer estimators ===========

#'@title PLS outer estimation
#'
#'@description
#'
#'Calculates a set of unstandardized outer weights.
#'
#'Mode A outer weights are correlations between indicators and composites.
#
#'Mode B outer weights are regression coefficients of composites on indicators. 
#'
#'For information about GSCA weights, see \link{GSCA}. 
#'
#'@inheritParams matrixpls-common
#'
#'@param ... All other arguments are ignored.
#'
#'@return A matrix of unscaled outer weights \code{W} with the same dimesions as \code{W.model}.
#'
#'@name outerEstimators
#'
#'@references
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial least squares.} Heidelberg: Physica-Verlag.
NULL

#'@describeIn outerEstimators Mode A outer estimation.
#'@export


outer.modeA <- function(S, W, E, W.model, ...){
  
  # Calculate the covariance matrix between indicators and composites
  W_new <- E %*% W %*% S
  
  # Set the non-existing weight relations to zero
  W_new[W.model == 0] <- 0
  
  return(W_new)
}

#'@describeIn outerEstimators Mode B outerestimation.
#'@export


outer.modeB <- function(S, W, E, W.model, ...){
  
  # Calculate the covariance matrix between indicators and composites
  IC <- E %*% W %*% S
  
  # Set up a weight pattern
  W_new <- ifelse(W.model==0,0,1)
  
  # Do the outer model regressions
  
  for(row in which(rowSums(W_new)>0, useNames = FALSE)){
    indicatorIndices <- W_new[row,]==1
    W_new[row,indicatorIndices] <- solve(S[indicatorIndices,indicatorIndices],IC[row,indicatorIndices])
  }
  
  return(W_new)
  
}

#'@describeIn outerEstimators outer estimation with generalized structured component analysis.
#@describeIn GSCA outer estimation with generalized structured component analysis.
#'@export


outer.gsca <- function(S, W, E, W.model, model, ...){
  
  nativeModel <- parseModelToNativeFormat(model)
  
  inner <- nativeModel$inner
  
  # Calculate the covariance matrix between indicators and composites and between composites
  IC <- W %*% S
  C <- W %*% S %*% t(W)
  
  # Estimate the reflective parts of the model
  
  reflective <- nativeModel$reflective
  
  for(row in which(rowSums(reflective!=0)>0)){
    independents <- which(reflective[row,] != 0)
    reflective[row,independents] <- solve(C[independents,independents],IC[independents, row])
  }  
  
  # Number of composites and indicators, and their sum
  P <- nrow(W.model)
  J <- ncol(W.model)
  JP <- J + P
  
  # E, and reflective form the A matrix of GSCA. 
  # Indicators first, then composites
  
  A <- rbind(reflective, E)
  V <- rbind(diag(J), W)
  
  # The following code is based on the ASGSCA package (licensed
  # under GPL-3). All matrices are transposed from the original
  # ASGSCA code
  
  # Step 2: Update W
  
  tr_w <- 0
  for(p in 1:P){
    t <- J + p
    windex_p <- which(W.model[p, ] != 0)
    m <- matrix(0, 1, JP)
    m[t] <- 1
    a <- A[, p]
    beta <- m - a
    H1 <- diag(P)
    H2 <- diag(JP)
    H1[p,p] <- 0
    H2[t,t] <- 0
    
    Delta <- A%*%H1%*%W - H2%*%V 
    Sp <- S[windex_p , windex_p]
    if (length(windex_p)!=0){        
      
      theta <- MASS::ginv(as.numeric(beta%*%t(beta))*S[windex_p,windex_p]) %*%
        t(beta %*% Delta %*% S[,windex_p])
      
      # Update the weights based on the estimated parameters and standardize
      W[p,windex_p] <- theta
      W <- scaleWeights(S, W)
      
    }
    
    # Proceed to next composite
  }
  
  
  return(W)
}
