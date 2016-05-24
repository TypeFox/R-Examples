H.fcnnls <-
function (x=x, y=y, verbose = FALSE, pseudo = FALSE, eps = 0)
{
# H.fcnnls : Function to calculate H matrix using Nonnegative Least Square method
# This function was obtained from NMF package (Gaujoux R., BMC Bioinformatics 2010,11:367) at
# https://github.com/renozao/NMF/blob/master/R/algorithms-snmf.R
# and was made minor bug fixes. 
# The original matlab code was proposed by M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450
# Given A and C this algorithm solves for the optimal 
# K in a least squares sense, using that
#      A = C*K 
# in the problem
#      min ||A-C*K||, s.t. K>=0, for given A and C.
# @param C the matrix of coefficients
# @param A the target matrix of observations
# @return [K, Pset]

# The function calculates H using estimated W and data
# x = W, y = dat
#--------------------------------------------------------------------------------
    if (any(dim(y) == 0L)) {
        stop("Empty target matrix 'y' [", paste(dim(y), collapse = " x "),
            "]")
    }
    if (any(dim(x) == 0L)) {
        stop("Empty regression variable matrix 'x' [", paste(dim(x),
            collapse = " x "), "]")
    }

    C <- x  
    A <- y  
    nObs = nrow(C)
    lVar = ncol(C)
    if (nrow(A) != nObs)
        stop("C and A have imcompatible sizes")
    pRHS = ncol(A)
    W = matrix(0, lVar, pRHS)
    iter = 0
    maxiter = 3 * lVar
    CtC = crossprod(C)        
    CtA = crossprod(C, A)     
    #K = .cssls(CtC, CtA, pseudo = pseudo)  # Original
K = adj.cssls(CtC, CtA)                
    Pset = K > 0

    K[!Pset] = 0
    D = K
    Fset = which(colSums(Pset) != lVar)
    oitr = 0
    while (length(Fset) > 0) {
        oitr = oitr + 1
        if (verbose && oitr > 5)
            cat(sprintf("%d ", oitr))
        K[, Fset] = adj.cssls(CtC, CtA[, Fset, drop = FALSE], Pset[,
        #    Fset, drop = FALSE], pseudo = pseudo)                    # Original
Fset, drop = FALSE])
        Hset = Fset[colSums(K[, Fset, drop = FALSE] < eps) >
            0]
        if (length(Hset) > 0) {
            nHset = length(Hset)
            alpha = matrix(0, lVar, nHset)
            while (nHset > 0 && (iter < maxiter)) {
                iter = iter + 1
                alpha[, 1:nHset] = Inf
                ij = which(Pset[, Hset, drop = FALSE] & (K[,
                  Hset, drop = FALSE] < eps), arr.ind = TRUE)
                i = ij[, 1]
                j = ij[, 2]
                if (length(i) == 0)
                  break
                hIdx = (j - 1) * lVar + i
                negIdx = (Hset[j] - 1) * lVar + i
                alpha[hIdx] = D[negIdx]/(D[negIdx] - K[negIdx])
                alpha.inf <- alpha[, 1:nHset, drop = FALSE]
                minIdx = max.col(-t(alpha.inf))
                alphaMin = alpha.inf[minIdx + (0:(nHset - 1) *
                  lVar)]
                alpha[, 1:nHset] = matrix(alphaMin, lVar, nHset,
                  byrow = TRUE)
                D[, Hset] = D[, Hset, drop = FALSE] - alpha[,
                  1:nHset, drop = FALSE] * (D[, Hset, drop = FALSE] -
                  K[, Hset, drop = FALSE])
                idx2zero = (Hset - 1) * lVar + minIdx
                D[idx2zero] = 0
                Pset[idx2zero] = FALSE
                K[, Hset] = adj.cssls(CtC, CtA[, Hset, drop = FALSE],
                #  Pset[, Hset, drop = FALSE], pseudo = pseudo)            # Original
  Pset[, Hset, drop = FALSE])
                Hset = which(colSums(K < eps) > 0)
                nHset = length(Hset)
            }
        }
        W[, Fset] = CtA[, Fset, drop = FALSE] - CtC %*% K[, Fset,
            drop = FALSE]
        Jset = which(colSums((ifelse(!(Pset[, Fset, drop = FALSE]),
            1, 0) * W[, Fset, drop = FALSE]) > eps) == 0)
        Fset = setdiff(Fset, Fset[Jset])
        if (length(Fset) > 0) {
            mxidx = max.col(t(ifelse(!Pset[, Fset, drop = FALSE],
                1, 0) * W[, Fset, drop = FALSE]))
            Pset[(Fset - 1) * lVar + mxidx] = TRUE
            D[, Fset] = K[, Fset, drop = FALSE]
        }
    }
    list(coef = K, Pset = Pset)
}
