GetRawCrCovFuncScal <- function(Ly, Lt = NULL, Ymu, Z,  Zmu ){
  # Sparse case if Ly and Lt are both lists
  if( is.list(Ly) && is.list(Lt) ){
    ulLt = unlist(Lt)
    if (length(Ymu) != length(unique(ulLt))){
      stop("Ymu and Lt are of the same size.")
    } else { 
      zstar = Z - Zmu;
      RCC <- list(tpairn = ulLt, 
                  rawCCov = rep(zstar, times = unlist( lapply(Ly,  length))) * 
                             (unlist(Ly) - approx(x= sort(unique(ulLt)), y = Ymu, xout = ulLt)$y ) )
      return(RCC)
    }
  # Dense case if Ly is a matrix and Lt is empty
  } else if ( is.matrix(Ly) && is.null(Lt)) {
    if( length(Z) != dim(Ly)[1] ) {
      stop("Ly and Z are not compatible (possibly different number of subjects).")
    }
    RCC <- list( tpairn = NULL,
                 rawCCov = cov( Ly, Z, use="pairwise.complete.obs" ))
  } else {
    stop("It appears you do no refine a valid cross-covariance type.")
  }
}

