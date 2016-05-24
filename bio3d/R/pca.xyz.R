"pca.xyz" <-
function(xyz, subset = rep(TRUE, nrow(as.matrix(xyz))), use.svd = FALSE,
         rm.gaps=FALSE, mass = NULL, ...) {
  ## Performs principal components analysis on the given "xyz" numeric data
  ## matrix and return the results as an object of class "pca.xyz"

  ## Log the call
  cl <- match.call()
  xyz <- as.xyz(xyz)

  if (any(!is.finite(xyz))) {
    ## Check for GAP positions in input
    if(rm.gaps) {
       gapC <- colSums(is.na(xyz)) == 0
      if (sum(gapC) > 3) {
        xyz <- xyz[,gapC]
        cat(paste("NOTE: Removing", sum(!gapC)/3, "gap positions with missing coordinate data\n",
            "     retaining", sum(gapC)/3, "non-gap positions for analysis.\n"))
      } else {
        stop("No non-gap containing positions (cols) available for analysis.")
      }
    } else {
       stop( paste("  Infinite or missing values in 'xyz' input.",
      "\t Likely solution is to remove gap positions (cols)",
      "\t or gap containing structures (rows) from input.", sep="\n") )
    }
  }

  dx <- dim(xyz)
  n <- dx[1]; p <- dx[2]
  if (!n || !p)
    stop("0 extent dimensions")
 
  # for mass-weighted PCA
  if(!is.null(mass)) {
     if(is.pdb(mass)) mass = aa2mass(mass)
     if(length(mass) != ncol(xyz)/3)
        stop("Input mass vector does not match xyz")
     q = t( t(xyz) * rep(sqrt(mass), each=3) )  # mass weighted xyz

     # re-do fitting: iteratively fit to the mean
     mean <- colMeans(q[subset, ])
     tolerance = 1.0 # convergence check
     maxiter = 10    # maximum number of iteration
     iter = 0
     repeat {
        q <- fit.xyz(mean, q, 1:ncol(q), 1:ncol(q), ...)
        mean.now <- colMeans(q[subset, ])
        mean.diff <- rmsd(mean, mean.now, 1:ncol(q), 1:ncol(q))
        mean = mean.now
        iter = iter + 1
        if(iter >= maxiter || mean.diff <= tolerance) break
     }
     if(mean.diff > tolerance) warning("Iteration stops before convergent")
     xyz <- q
  }
 
#  mean <- apply(xyz[subset,],2,mean) ## mean structure
  mean <- colMeans(xyz[subset,]) ## Faster
  n <- sum(subset) 
  
  # Check number of columns
  if(p > 3000 && n <= 0.4*p && !use.svd) {
     cat("NOTE: In input xyz (MxN),  N > 3000 and M < N\n",
         "     Singular Value Decomposition (SVD) approach is faster\n",
         "     and is recommended (set 'use.svd = TRUE')\n\n", sep=" ")
     flush(stdout())
  }
     
  if(!use.svd) {   
     S    <- var(xyz[subset,])          ## coverance matrix
   
     ## eigenvectors ("U") & eigenvalues ("L"): [ U'SU=L ]
     prj  <- eigen(S, symmetric = TRUE)
     L <- prj$values
     U <- prj$vectors
  } else {
     if(n < p)
        warning(paste("In input xyz (MxN), M < N:\n",
           "   Only",n,"eigenvalues and eigenvectors are returned!\n\n"))

     ## S = Q'Q, Q = UDV'
     Q <- t(t(xyz[subset,]) - mean) / sqrt(n-1)
     prj <- svd(Q)
     L <- prj$d^2
     U <- prj$v
  }

  ## fix negative eigenvalues
  ## (these are very small numbers and should be zero)
  L[L<0]<-0
  sdev <- sqrt(L)

  ## scores of "xyz" on the pc's [ z=U'[x-x.mean] ]
  z <- sweep(xyz,2,mean) %*% (U)

  ## atom-wise loadings (norm of xyz eigenvectors)
  ## Skip the calculation if the input is not xyz coordinates,
  ## e.g. for PCA over correlaiton matrices (see pca.array()).
  if(ncol(U) %% 3 == 0) { 
     au <- apply(U, 2, function(x) {
       sqrt(colSums(matrix(x^2, nrow=3))) })
  } else {
     au <- NULL
  }
  
  class(U)="pca.loadings"

  if(!is.null(mass)) {
     mean = mean / sqrt(rep(mass, each=3))
     out <- list(L=L, U=U, z=z, au=au,
              sdev=sdev, mean=mean, mass=mass, call=cl)
  }
  else
     out <- list(L=L, U=U, z=z, au=au,
              sdev=sdev, mean=mean, call=cl)

  class(out)="pca"; out
}
