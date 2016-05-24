################################################################################
#  
#  Extension for gam package allowing a projection Slepian basis to be used;
#  part of the slp package
#  (C) wburr, July 2014
#
#  Licensed under GPL-2
#
################################################################################

slp <- function(x, W = NA, K = NA, deltat = 1, naive = FALSE, 
                intercept = FALSE, customSVD = TRUE, forceC = FALSE) {

  # logical checks (note: W is assumed to be in the same units as deltat, i.e., days)
  stopifnot(is.numeric(x), (is.na(W) | is.numeric(W) & W > 0 & W < 1/(2*deltat)), deltat %in% c(1, 6),
            (is.numeric(K) & K <= length(x) | is.na(K)), is.logical(customSVD))

  namesx <- names(x)
  x <- as.vector(x)

  # x should be a contiguous time array; if not, convert to one, back-convert at the end
  if(is.na(max(x[-1L] - x[-length(x)])) | max(x[-1L] - x[-length(x)]) > deltat) {
    warning("slp: Input time array is not contiguous.")
    minT <- min(x, na.rm = TRUE); maxT <- max(x, na.rm = TRUE)
    stopifnot(is.numeric(minT), is.numeric(maxT))
    wx <- seq(minT, maxT, deltat)

    if(deltat == 6) {
      wx <- seq(minT, maxT, deltat)
      if(wx[length(wx)] != maxT) { stop("slp: Input time array not properly time aligned to 6-day samples per deltat = 6.") }
    }

          cat(" Input Time Array: \n")
    cat(paste(" *        samples: ", length(x), "\n", sep = ""))
    cat(paste(" *      max - min: ", maxT - minT + 1, "\n", sep = ""))
    cat(paste(" *         deltat: ", deltat, "\n", sep = ""))
  } else {
    wx <- x 
  }

  N <- length(wx)


  if(is.na(W) & is.na(K)) { stop("Must set one of K or W for family selection.") }
  if(!is.na(K) & !is.na(W)) {
      if(K > ceiling(2 * N * W)) { stop("Using more than 2NW basis vectors is not recommended (or supported).") }
  } else {
      if(is.na(K)) {
          K <- as.integer(round(2 * N * W))
      } else {   # Case of is.na(W), K is set
        if(floor(K) != K) {  # K doesn't have to be integer _class_, but it needs to be an integer
          K <- floor(K)
          warning(paste("slp: K choice not integer. Truncated to K = ", K, sep = ""))
        }
        W <- as.numeric((K)) / as.numeric((2 * length(wx)))
      }
  }

  # Logical check: if N, K and W match one of the saved objects, simply
  # return that object; otherwise, run through the generation process
  Wn <- round(W * 365.2425)    # convert to df/year
  if(checkSaved(N, Wn, K) & !forceC) {

    data(list = as.character(paste0("basis_N_", N, "_W_", Wn, "_K_", K)), 
         envir = environment())

    if(!intercept) { basis <- basis[, -1] }
    
    # need to convert back to original time array -- the NAs in the original
    if(any(is.na(x))) {
      basis[which(!(wx %in% x[!is.na(x)])), ] <- rep(NA, ncol(basis))
    }
  } else {

      # start by generating baseline Slepian basis vectors
      v <- .dpss(n = N, nw = N * W, k = K)
      if(naive) {        # Case of non-mean-adjusted Slepians, SLP
        basis <- v
      } else {           # Case of mean-adjusted Slepians, SLP2 or SLP3
    
        # Equations follow from Thomson (2001), see help file for full citation
        alpha <- t(rep(1, N)) %*% v %*% t(v) %*% rep(1, N) 
        R <- rep(1/sqrt(alpha), N)         
        U <- t(v) %*% R                    

        #  Two options: intercept = TRUE / FALSE
        sRaw <- v %*% (diag(K) - U %*% t(U)) %*% t(v)

        if(customSVD) {
          basis <- .slpsvd(A = sRaw, N = N, K = K)$u      # optimized LAPACK SVD for our edge case
        } else { 
          basis <- svd(sRaw, nu = K, nv = 0)$u            # generic R -> C -> LAPACK SVD
        }

        basis <- basis[, -K]   # mean adjustment implies a mean-passing subspace of dim K-1,
                               # conveniently arranged so that basis vector K is out-of-band
        if(intercept) { basis <- cbind(R, basis) }

        dimnames(basis) <- list(names(wx), 1L:ncol(basis))
     
        # need to convert back to original time array -- the NAs in the original
        if(any(is.na(x))) {
          basis[which(!(wx %in% x[!is.na(x)])), ] <- rep(NA, ncol(basis))
        }
  
        a <- list(K = K, W = W, N = N, naive = naive)
        attributes(basis) <- c(attributes(basis), a)
        class(basis) <- c("slp", "basis", "matrix")
      }
  }
  basis
}

