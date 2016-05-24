CoxFindU3 <- function(BetaMLE, StepSize, Hfun, Efun, y, d, Z, level=3.84) 
{
    MaxValue <- CoxEL(y=y, d=d, Z=Z, beta=BetaMLE, lam=0, fun=Hfun)$logEmpLik
    target <- (MaxValue - level/2)
    print(target)
    stepsize <- StepSize
    b1vec <- BetaMLE[1] + stepsize[1] * (1:10 - 5.5)
    b2vec <- BetaMLE[2] + stepsize[2] * (1:10 - 5.5)
    lamvec <- stepsize[3] * (1:10 - 5.5)

    temp <- matrix(NA, nrow = 5, ncol = 1000)
    temp[1, ] <- rep(b1vec, each = 100)
    temp[2, ] <- rep(rep(b2vec, each = 10), 10)
    temp[3, ] <- rep(lamvec, 100)
    for (i in 1:10) for (j in 1:10) for (k in 1:10) {
        ind <- (i - 1) * 100 + (j - 1) * 10 + k
        mtemp <- CoxEL(y=y, d=d, Z=Z, beta=c(b1vec[i], b2vec[j]), lam=lamvec[k], fun=Hfun)
        temp[4, ind] <- mtemp$mu
        temp[5, ind] <- mtemp$logEmpLik
    }

    subsetInd <- ( temp[5, ] >= target )
    ParaTrials <- temp[, subsetInd]
    TrialValues <- Efun(beta=cbind(ParaTrials[1,], ParaTrials[2,]), theta=ParaTrials[4,])
    maxInd <- which.max(TrialValues)
    EfunU <- TrialValues[maxInd]
    maxPara <- ParaTrials[, maxInd]
    print(c(EfunU, maxPara))


    for (N in 1:30) {
        if (sum(subsetInd) < 500) 
            stepsize <- stepsize/4
        b1vec <- maxPara[1] + stepsize[1] * (1:10 - 5.5)
        b2vec <- maxPara[2] + stepsize[2] * (1:10 - 5.5)
        lamvec <- maxPara[3] + stepsize[3] * (1:10 - 5.5)
        temp[1, ] <- rep(b1vec, each = 100)
        temp[2, ] <- rep(rep(b2vec, each = 10), 10)
        temp[3, ] <- rep(lamvec, 100)
        for (i in 1:10) for (j in 1:10) for (k in 1:10) {
            ind <- (i - 1) * 100 + (j - 1) * 10 + k
            mtemp <- CoxEL(y=y, d=d, Z=Z, beta=c(b1vec[i], b2vec[j]), lam=lamvec[k], fun=Hfun)
            temp[4, ind] <- mtemp$mu
            temp[5, ind] <- mtemp$logEmpLik
        }

        subsetInd <- (temp[5, ] >= target)
        ParaTrials <- temp[, subsetInd]
        TrialValues <- Efun(beta=cbind(ParaTrials[1,],ParaTrials[2,]), theta=ParaTrials[4,])
        maxInd <- which.max(TrialValues)
        EfunU <- TrialValues[maxInd]
        maxPara <- ParaTrials[, maxInd]
        print(c(EfunU, maxPara))
    }
    list(Upper = EfunU, maxParameterNloglik = maxPara)
}