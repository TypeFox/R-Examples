nortestARMA <- function(residuals, sig2hat, d = 2) {

    if (d > 10) stop("'d' should be an integer smaller than 10.")
    nT <- length(residuals)
    RK <- RKsdest <- RKsdmuest <- rep(0, 10)
    Kopt <- Koptsdest <- Koptsdmuest <- 0
    st <- stsdest <- stsdmuest <- 0
    pval <- pvalsdest <- pvalsdmuest <- 0

    out <- .C("calcstat", residuals = as.double(residuals), sigmahat = as.double(sqrt(sig2hat)), nT = as.integer(nT), d = as.integer(d), RK = as.double(RK),
              RKsdest = as.double(RKsdest), RKsdmuest = as.double(RKsdmuest), Kopt = as.integer(Kopt), Koptsdest = as.integer(Koptsdest), Koptsdmuest = as.integer(Koptsdmuest), st = as.double(st), stsdest = as.double(stsdest), stsdmuest = as.double(stsdmuest), pval = as.double(pval), pvalsdest = as.double(pvalsdest), pvalsdmuest = as.double(pvalsdmuest),
                    PACKAGE = "nortestARMA")

    names(out$RKsdest) <- 1:length(out$RKsdest)
    names(out$RKsdmuest) <- 1:length(out$RKsdmuest)
    
    return(list(RKsdest = out$RKsdest, Koptsdest = out$Koptsdest, pvalsdest = out$pvalsdest, RKsdmuest = out$RKsdmuest, Koptsdmuest = out$Koptsdmuest, pvalsdmuest = out$pvalsdmuest))
 
}

