########################################################################
#
#   Extensions to mgcv package allowing a projection Slepian basis
#   to be used; part of the slp package
# 
#   (C) wburr, July 2014
#   Licensed under GPL-2
#
########################################################################


########################################################################
#
#   Smooth Constructor
#
#   Forms model matrix from basis vectors consisting of orthogonal
#   Slepian sequences (DPSSs). 
#
########################################################################
smooth.construct.slp.smooth.spec <- function(object, data, knots) {

    # p.order is inapplicable; we are not using a polynomial
    if (!is.na(object$p.order[1])) warning("Specifying 'p' is meaningless for slp smoothers.")
  
    # pull parameters out of 'xt'
    ext <- object[['xt']]
  
    # W, if specified
    if(!is.null(ext[['W']])) { 
        W <- ext[['W']] 
        stopifnot(is.numeric(W), W > 0, W < 0.5)
    } else { W <- NA }
  
    # K, if specified
    if(!is.null(ext[['K']])) {
        K <- ext[['K']]
        if(floor(K) != K) {  # K doesn't have to be integer _class_, but it needs to be an integer
            K <- floor(K)
            warning(paste("slp: K choice not integer. Truncated to K = ", K, sep = ""))
        }
    } else {
        K <- NA
    }
  
    # check for bs.dim setting -- if left blank, mgcv defaults to -1
    if(object[['bs.dim']] != -1) {
        warning("slp: Set dimension via xt() as 'K'. Parameter ignored.")
    }
  
    # deltat, if specified (for 6-day-sampled data)
    if(!is.null(ext[['deltat']])) {
        deltat <- ext[['deltat']]
        stopifnot(round(deltat) == deltat, deltat %in% c(1, 6))
    } else {   # defaults to 1
        deltat <- 1
    }
  
    # x array
    x <- data[[object$term]]
    if(is.na(max(x[-1L] - x[-length(x)])) | max(x[-1L] - x[-length(x)]) > deltat) {
        warning("slp: specified input time array is not contiguous (missing values). Correcting.")
        contiguous <- FALSE
  
        if(!is.null(ext[['mask']])) {
            mask <- ext[['mask']]
            isMask <- TRUE
  
            # sanity check 1: is mask TRUE/FALSE?
            if(length(c(which(mask == TRUE), which(mask == FALSE))) != length(mask) ) {
                stop("'mask' must be populated with TRUE/FALSE elements.")
            }
  
            # sanity check 2: mask specified, does it match up with object$term?
            if(length(which(mask == TRUE)) != length(x)) {
                cat(paste0("\n", str(x), " (Data)", "\n"))
                cat(paste0(str(mask), " (Mask)", "\n"))
                cat(paste0(length(x), " good data points vs ", length(which(mask == TRUE)), " true elements \n"))
                stop("Mask must correspond to missing data.")
            }
            wx <- seq(from = 1, length.out = length(mask), by = deltat)
  
        } else {
            warning("slp: mask variable not specified; interpolating time.")
            isMask <- FALSE
  
            minT <- min(x, na.rm = TRUE); maxT <- max(x, na.rm = TRUE)
            stopifnot(is.numeric(minT), is.numeric(maxT))
            wx <- seq(minT, maxT, deltat)
  
            if(deltat == 6) {
                wx <- seq(minT, maxT, deltat)
                if(wx[length(wx)] != maxT) { 
                    stop("slp: Input time array not properly time aligned to 6-day samples per deltat = 6.") 
                }
            }
                  cat("slp: Input Time Array: \n")
            cat(paste(" *            samples: ", length(x), "\n", sep = ""))
            cat(paste(" *          max - min: ", maxT - minT + 1, "\n", sep = ""))
            cat(paste(" *             deltat: ", deltat, "\n", sep = ""))
        }
      # end of "some missing data, figure out time index array"
    } else {
        contiguous <- TRUE
        isMask <- FALSE
        wx <- x 
    }
    
    # N, if specified
    if(!is.null(ext[['N']])) {   
        N <- ext[['N']]
        stopifnot(round(N) == N, N > 1)
    } else {
        N <- length(wx)
    }
  
    if(is.na(W) & is.na(K)) { stop("Must set one of K or W for family selection.") }
    if(!is.na(K) & !is.na(W)) {
        if(K > ceiling(2 * N * W)) { stop("Using more than 2NW basis vectors is not recommended (or supported).") }
    } else {
        if(is.na(K)) {
            K <- as.integer(round(2 * N * W))
        } else {   # Case of is.na(W)
            W <- as.numeric((K)) / as.numeric((2 * length(wx)))
        }
    }
  
    # At this point: K, W and N should be specified properly
    stopifnot(round(K) == K, is.numeric(W), round(N) == N, K > 1, W < 0.5, W > 0)
  
    # Other, non-numeric, parameters:
  
    # naive (SLP vs SLP2 or SLP3)
    if(!is.null(ext[['naive']])) {
        naive <- ext[['naive']]
    } else {
        naive <- FALSE
    }
  
    # intercept, choosing SLP2 vs SLP3
    if(!is.null(ext[['intercept']])) {
        intercept <- ext[['intercept']]
    } else {
        intercept <- FALSE
    }
  
    # use default SVD (slow) or customSVD (faster); should return same result regardless of choice
    if(!is.null(ext[['customSVD']])) {
        customSVD <- ext[['customSVD']]
    } else {
        customSVD <- TRUE
    }

    # force computation of basis vectors, or read from disk (if available)?
    if(!is.null(ext[['forceC']])) {
        forceC <- ext[['forceC']]
    } else {
        forceC <- FALSE
    }

    ################################################################################   
    # 
    #  All parameters set: 3 numeric (N, W, K), 3 logical (naive, intercept, customSVD),
    #  possibly mask, possibly forceC
    #  
    #  From here the code is very similar to slp-gam.R, with the only change being
    #  in the careful examination of a 'mask' variable _if_ wx != data[['term']]
    #
    ################################################################################
  
    #
    #  Notes:
    #  * the passed in 'object' is the return
    #  * the parameter 'C' in object determines the centering constraints
    #    (intercept impacts this)
    #    * the Slepian functions are orthonormal, and the SVD'd SLP2 and SLP3 bases
    #      are also in good condition, and should not be messed with
    #    * by default, mgcv will normalize them, and form a QR decomposition
    #    * this is not desirable in this case, so we set C = zero-row matrix
    #
    object[['C']] <- matrix(data=NA, nrow = 0, ncol = K)

    Wn <- round(W * 365.2425)
    if(checkSaved(N, Wn, K) & !forceC) {  # this case is here to load the basis set from hard disk (saved)

        data(list = as.character(paste0("basis_N_", N, "_W_", Wn, "_K_", K)),
             envir = environment())
 
        if(!intercept) { basis <- basis[, -1] }
  
        basisFull <- basis
  
        # need to convert back to original time array -- the NAs in the original
        if(!contiguous) {
            if(isMask) {
                basis <- basis[which(mask == TRUE), ]
            } else {
                basis <- basis[which(wx %in% x), ] 
            }
        }
    } else {
  
        # start by generating baseline Slepian basis vectors
        v <- .dpss(n = N, nw = N * W, k = K)
  
        if(naive) {        # Case of non-mean-adjusted Slepians, SLP
            basis <- v
            basisFull <- basis
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
            basisFull <- basis 
    
            # need to convert back to original time array -- the NAs in the original
            if(!contiguous) {
                if(isMask) {
                    basis <- basis[which(mask == TRUE), ]
                } else {
                    basis <- basis[which(wx %in% x), ] 
                }
            }
        } # end of SLP2/SLP3

        a <- list(K = K, W = W, N = N, naive = naive)
        attributes(basis) <- c(attributes(basis), a)
        class(basis) <- c("slp", "basis", "matrix")
    } # end of basis vector generation
  
    object[['size']] <- c(length(basis[, 1]), length(basis[1, ]))
    object[['X']] <- basis
    object[['rank']] <- length(basis[1, ])
    object[['null.space.dim']] <- object[['size']][1] - object[['size']][2]
  
    # store 'slp' specific stuff
    if(isMask) { object[['mask']] <- mask }
    object[['K']] <- K
    object[['W']] <- W
    object[['N']] <- N
  
    # slp does not support penalization, so ...
    object[['fixed']] <- TRUE
    object[['bs.dim']] <- object[['size']][2]
    object[['C']] <- object[['C']][, 1:object[['bs.dim']]]   # correct for number of vectors in final basis
    object[['fullBasis']] <- basisFull
    object[['contiguous']] <- contiguous
    object[['wx']] <- wx

    class(object)<-"slp.smooth"  # Give object a class
    object
}


########################################################################
#
#   Predictor
#
#   Forms prediction (for summary, predict, model.matrix, etc.) from 
#   model matrix using provided (sub-set) X. 
#
########################################################################
Predict.matrix.slp.smooth<-function(object,data)
{ 
    # Peel the required basis out of the object; for some reason, the computed
    # basis vectors aren't passed in, hence the necessity of saving the 'fullBasis'
    # in the constructor above -- they're computationally intensive!

    px <- data[[object[['term']]]]

    # full basis object is saved as part of the object (large, but it saves computation)
    basis <- object[['fullBasis']]
    contiguous <- object[['contiguous']]
    isMask <- !is.null(object[['mask']])
    wx <- object[['wx']]

    if(!contiguous) {
        if(isMask) {
            X <- basis[which(object[['mask']] == TRUE), ]
        } else {
            X <- basis[which(wx %in% px), ] 
        }
    } else {
        X <- basis
    }

    X
}

