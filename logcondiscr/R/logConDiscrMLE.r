logConDiscrMLE <- function(x, w = NA, psi_o = NA, prec = 1e-5, output = TRUE){

########################################################################
#
# Compute discrete, concave function psi on [X(1), X(n)] with knots
# only in {X(1), X(2),...,X(m)} such that the log-likelihood function
#     L(psi) <- sum_{i=1}^{m} w(i)*psi(X(i)) 
#              - sum_{j=X(1)}^{X(m)} exp(psi(j))
# is maximal.
#
# INPUT
# - x      : vector of strictly increasing observations
# - w      : vector of weights with entries W_k:= #{i : XRawData_i == X_k}/n
# - psi_o  : concave initialization vector
#
# OUTPUT
# - psi    : vector with entries psi_i := psi(X_i)
# - L      : value of log-likelihood function at maximum
# - isKnot : 0-1-vector where isKnot_k = 1 if psi has a knot at X_k
# - psiSupp: psi interpolated on {X_1, ..., X_m] 
#
# Notwendige m-Files:
# dMLE.m (deshalb auch J00.m, J10.m, J11.m, J20.m)
#
# Matlab version 20.11.07, Kathrin Weyermann
# Ported to R by Kaspar Rufibach, October 2010
########################################################################

# prepare data
if (identical(w, NA)){
    X <- unique(sort(x))
    W <- as.vector(table(x) / length(x))
    }
    
if (identical(w, NA) == FALSE){
    w <- as.numeric(w)
    W <- w[order(x)]
    X <- sort(x)
    W <- w / sum(w)
    }    

# initalize
iterActive1 <- 0
dX <- diff(X)
m <- length(X)
isKnot <- rep(0, m)
if (identical(psi_o, NA)){psi_o <- LocalNormalize(rep(1, m), dX)} else {
    psi_o <- LocalNormalize(psi_o, dX)  
    conc <- LocalConcavity(psi_o, dX)
    isKnot <- as.numeric(abs(conc) > prec)
    }
isKnot[1] <- 1 
isKnot[length(isKnot)] <- 1

# active set algorithm
res <- LocalMLE(X, W, isKnot, psi_o, prec)
psi <- res$psi
L <- res$L
conc <- res$conc
D <- res$D

# basic procedure 2
while (((max(D) > prec * mean(abs(D)))) & (iterActive1 < 500)){
    iterActive1 <- iterActive1 + 1
    isKnot_old <- isKnot
    #Knoten hinzufuegen:
    k <- which.max(D)
    isKnot[k] <- 1    
    res <- LocalMLE(X, W, isKnot, psi, prec)
    psi_new <- res$psi
    L <- res$L
    conc_new <- res$conc
    D <- res$D    
    
    # basic procedure 1
    iterActive2 <- 0
    while (max(conc_new) > prec * max(abs(conc_new))){

        # remove non-concave knots
        JJ <- which(conc_new > 0)  
        conc_diff <- conc[JJ] - conc_new[JJ]
        
        # avoid division by 0
        conc_diff[conc_diff == 0] <- prec ^ 2    
        tmp <- conc[JJ] / conc_diff
        lambda <- min(tmp)
        KK <- which(tmp == lambda)
        isKnot[JJ[KK]] <- 0
        psi <- (1 - lambda) * psi + lambda * psi_new
        conc <- pmin(LocalConcavity(psi, dX), 0)
        res <- LocalMLE(X, W, isKnot, psi, prec)
        psi_new <- res$psi
        L <- res$L
        conc_new <- res$conc
        D <- res$D
        iterActive2 <- iterActive2 + 1
        
        if (identical(output, TRUE)){cat(paste("Outer active set iteration:", format(iterActive1, width = 3), " / inner iteration:", format(iterActive2, width = 2), "\n", sep = ""))}
    } # end while
    psi <- psi_new
    conc <- conc_new
    
    # avoid endless loop if set of knots does not change
    if (identical(isKnot, isKnot_old)){D <- rep(0, length(D))}

    if (identical(output, TRUE)){cat(paste("Outer active set iteration:", format(iterActive1, width = log10(iterActive1) + 2), "\n", sep = ""))}
    } # end while
    
# interpolate estimate on {X_1, ..., X_m}
II <- NULL
XSupp <- X[1]:X[length(X)]   
for (i in 1:length(X)){II <- c(II, which(XSupp == X[i]))}
isObs <- rep(0, length(XSupp))
isObs[II] <- 1
psiSupp <- LocalExtend(XSupp, isObs, X, psi)

# collect and output results
res <- list("x" = X, "w" = W, "psi" = psi, "L" = L, "isKnot" = isKnot, "xSupp" = XSupp, "psiSupp" = psiSupp)
return(res)
}
