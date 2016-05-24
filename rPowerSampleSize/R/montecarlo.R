montecarlo <- function(method, M = 100000, nE, r, m, nCovernE = 1, muC, muE, d = rep(0.0, m), SigmaE, SigmaC, alpha = 0.05, q = 1, nbcores = parallel::detectCores() - 1, alternative = "greater", orig.Hochberg = FALSE) {

    if (!(method %in% c("Bonferroni", "Hochberg", "Holm"))) stop("'method' argument misspecified.")
    if (orig.Hochberg && (q != 1)) stop("'orig.Hochberg' should be set to TRUE only when q=1.") 

    matrix.type <- matrix.type.compute(SigmaE, SigmaC, display.type = FALSE)
    if ((matrix.type == 3) || (matrix.type == 4)) rho <- SigmaE[1, 2] / SigmaE[1, 1] else rho <- NULL
    equalSigmas <- isTRUE(all.equal(SigmaE, SigmaC))
    
    p <- m
    nC <- nCovernE * nE

#    D1qm <- .D1(q, m)
#    alphai <- rep(NA, m)
#    for (i in 1:q) alphai[i] <- q / m
#    if ((q + 1) <= m) for (i in (q + 1):m) alphai[i] <- q / (m + q - i)

    myfunc <- function(M, method, d, matrix.type, equalSigmas, alpha, q, rho, alternative, nE, nC, muE, muC, SigmaE, SigmaC) {
        z <- rep(NA, M)
        for (b in 1:M) {
            XE <- mvtnorm::rmvnorm(nE, mean = muE, sigma = SigmaE)
            XC <- mvtnorm::rmvnorm(nC, mean = muC, sigma = SigmaC)
            res <- indiv.analysis(method = method, XE = XE, XC = XC, d = d, matrix.type = matrix.type, equalSigmas = equalSigmas, alpha = alpha, q = q, rho = rho, alternative = alternative, orig.Hochberg = orig.Hochberg)
                                        #        pvals <- sort(res$pvals)
                                        #        nb <- 0
                                        #        for (j in m:1) {
                                        #            if (pvals[j] < alpha * alphai[j] / D1qm) {nb <- j; break}
                                        #        }
                                        #        z[b] <- nb            
            z[b] <- sum(res$AdjPval < alpha)
        }
        return(z)
#        return(sum(z >= r))
    }

    if (nbcores > 1) { # We start the cluster
        cl <- parallel::makeCluster(getOption("cl.cores", nbcores))	
        out <- parallel::clusterCall(cl, myfunc, M = round(M / nbcores), method = method, d = d, matrix.type = matrix.type, equalSigmas = equalSigmas, alpha = alpha, q = q, rho = rho, alternative = alternative, nE = nE, nC = nC, muE = muE, muC = muC, SigmaE = SigmaE, SigmaC = SigmaC) # M/nbcores iterations are performed on each core          
      # We stop the cluster
        parallel::stopCluster(cl)
        vecz <- unlist(out)
        res <- sum(vecz >= r) / (nbcores * round(M/nbcores))
    } else {
        vecz <- myfunc(M, method = method, d = d, matrix.type = matrix.type, equalSigmas = equalSigmas, alpha = alpha, q = q, rho = rho, alternative = alternative, nE = nE, nC = nC, muE = muE, muC = muC, SigmaE = SigmaE, SigmaC = SigmaC)
        res <- sum(vecz >= r) / M
    }

    
    if (method == "Bonferroni") result <- list(rpowBonf = res, vecr = vecz)
    if (method == "Hochberg") result <- list(rpowHoch = res, vecr = vecz)
    if (method == "Holm") result <- list(rpowHolm = res, vecr = vecz)

    class(result) <- "rPower"

    return(invisible(result))
}

