FitARz <-
function (z, p, demean = TRUE, MeanMLEQ = FALSE, lag.max = "default") 
{
    stopifnot(length(z) > 0, length(z) > max(p), length(p) > 
        0)
    is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    stopifnot(is.wholenumber(p), p>0)
    ztsp <- tsp(z)
    if (lag.max == "default") 
        MaxLag <- min(300, ceiling(length(z)/5))
    else    MaxLag = lag.max
    MaxIter <- 10
    n <- length(z)
    pvec <- sort(p)
    pvec <- pvec[pvec > 0]
    if (length(pvec) == 0) 
        pvec <- 0
    if (length(p) == 1 && pvec != 0) 
        pvec <- 1:p
    PMAX <- max(pvec)
    SubQ <- length(pvec) < PMAX
    indMeanQ <- demean || MeanMLEQ
    if (indMeanQ) 
        mz <- mean(z)
    else mz <- 0
    y <- z - mz
    ans <- GetFitARz(y, pvec)
    LL <- ans$loglikelihood
    etol <- 1
    mu <- iter <- 0
    if (MeanMLEQ && PMAX != 0) 
        while (etol > 1e-06 && iter < MaxIter) {
            LLPrev <- LL
            iter <- iter + 1
            mu <- GetARMeanMLE(y, ans$phiHat)
            ans <- GetFitAR(y - mu, pvec)
            LL <- ans$loglikelihood
            etol <- abs(LL - LLPrev)/LLPrev
            if (ans$convergence != 0) 
                stop("GetARFit returned convergence = ", ans$convergence)
        }
    muHat <- mu + mz
    zetaHat <- ans$zetaHat
    phiHat <- ans$phiHat
    if (PMAX != 0) 
        res <- BackcastResidualsAR(y, phiHat, Q = 100, demean = FALSE)
    else res <- y
    fits <- y - res
    sigsq <- sum(res^2)/n
    racf <- (acf(res, plot = FALSE, lag.max = MaxLag)$acf)[-1]
    if (SubQ) {
        varNames <- paste("zeta(", pvec, ")", sep = "")
        covHat <- solve(InformationMatrixARz(zetaHat, pvec))/n
        dimnames(covHat) <- list(varNames, varNames)
        sdRacf <- sqrt(diag(VarianceRacfARz(zetaHat, pvec, MaxLag, 
            n)))
    }
    else {
        if (PMAX > 0) {
            varNames <- paste("phi(", 1:PMAX, ")", sep = "")
            covHat <- SiddiquiMatrix(phiHat)/n
            dimnames(covHat) <- list(varNames, varNames)
            sdRacf <- sqrt(diag(VarianceRacfAR(phiHat, MaxLag, 
                n)))
        }
        else {
            varNames <- character(0)
            covHat <- numeric(0)
            sdRacf <- rep(1/sqrt(n), MaxLag)
        }
    }
    RacfMatrix <- matrix(c(racf, sdRacf), ncol = 2)
    dimnames(RacfMatrix) <- list(1:MaxLag, c("ra", "Sd(ra)"))
    LBQ <- LjungBoxTest(res, lag.max = MaxLag, k = length(zetaHat))
    if (SubQ) {
        m <- length(pvec)
        if (m < 13) {
            pVEC <- deparse(as.numeric(pvec), width.cutoff = 180)
            pVEC <- substr(pVEC, 2, nchar(pVEC))
            }
        else {
            pVECa <- deparse(as.numeric(pvec[1:4]), width.cutoff = 180)
            pVECa <- substr(pVECa, 2, nchar(pVECa)-1)
            pVECb <- deparse(as.numeric(pvec[(m-2):m]), width.cutoff = 180)
            pVECb <- substr(pVECb, 3, nchar(pVECb)) 
            pVEC <- paste(pVECa, ",...,", pVECb, ", m=",m)
         } 
        ModelTitle <- paste("ARz", pVEC, sep = "")
        ModelTitle <- gsub(" ", "", ModelTitle)
    }
    else ModelTitle <- paste("AR(", p, ")", sep = "")
    ans <- list(loglikelihood = ans$loglikelihood, phiHat = phiHat, 
        sigsqHat = sigsq, muHat = muHat, covHat = covHat, zetaHat = zetaHat, 
        RacfMatrix = RacfMatrix, LjungBoxQ = LBQ, res = res, 
        fits = fits + mz, SubsetQ = SubQ, pvec = pvec, demean = demean, 
        FitMethod = "MLE", iterationCount = iter, convergence = ans$convergence, 
        MeanMLE = MeanMLEQ, tsp = ztsp, call = match.call(), 
        ARModel = "ARz", DataTitle = attr(z, "title"), ModelTitle = ModelTitle, 
        z = z)
    class(ans) <- "FitAR"
    ans
}

