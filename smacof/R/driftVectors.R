## embedding skew-symmetries as drift vectors

driftVectors <- function(data, type = c("ratio", "interval", "ordinal","mspline"), 
                         weightmat = NULL, init = "torgerson", ties = "primary",  verbose = FALSE, 
                         relax = FALSE, modulus = 1, itmax = 1000, eps = 1e-6, 
                         spline.degree = 2, spline.intKnots = 2)  {
  
## data ... asymmetric dissimilarity matrix
  
  symres <- symdecomp(data)   ## decompose data into symmetric and skew-symmetric portion
  M <- symres$M               ## symmetric matrix
  N <- symres$N               ## skew-symmetric matrix
  n <- nrow(M)                ## number of objects

  fitsym <- mds(M, ndim = 2, type = type, weightmat = weightmat, init = init, ties = ties, verbose = verbose,
                relax = relax, modulus = modulus, itmax = itmax, eps = eps, 
                spline.degree = spline.degree, spline.intKnots = spline.intKnots)
  x <- fitsym$conf
  ind <- expand.grid(1:n, 1:n)[,2:1]
  indmat <- ind[-which(ind[,1] == ind[,2]),]

  Amat <- t(apply(indmat, 1, function(ij) x[ij[2],] - x[ij[1],]))  ## conf difference
  Bmat <- t(apply(Amat, 1, function(ab) ab/sqrt((t(ab) %*% ab))))  ## unit length
  diag(N) <- NA
  nij <- as.vector(na.omit(c(N)))
  Cmat <- Bmat * nij
  rownames(Cmat) <- paste0("c", apply(indmat, 1, paste0, collapse = ""))
  Dmat <- apply(Cmat, 2, function(cc) tapply(cc, indmat[,1], mean))
  Dlength <- apply(Dmat, 1, function(dd) sqrt(sum(dd^2)))                     ## vector length
  names(Dlength) <- rownames(x)
  Dlength <- Dlength*n                                          ## due to normalization
  u <- c(0,1)
  alpha <- apply(Dmat, 1, function(di) acos((t(di) %*% u)/sqrt(t(di) %*% di)))  ## angle (radians)
  driftcoor <- cbind(x[,1] + cos(alpha) * Dlength, x[,2] - sin(alpha) * Dlength) ## drift coordinates
  result <- list(fitsym = fitsym, sym = M, skewsym = N, driftcoor = driftcoor, nobj = n, stress = fitsym$stress, 
                 niter = fitsym$niter, call = match.call())
  class(result) <- "driftvec"
  return(result)
}

