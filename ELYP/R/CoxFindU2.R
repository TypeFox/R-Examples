CoxFindU2 <- function(BetaMLE, StepSize, Hfun, Efun, y, d, Z, level=3.84) 
{
    MaxValue <- CoxEL(y=y, d=d, Z=Z, beta=BetaMLE, lam=0, fun=Hfun)$logEmpLik
    target <- MaxValue - level/2
    print(target)
    stepsize <- StepSize
    b1vec <- BetaMLE + stepsize[1]*(1:10 - 5.5)
    lamvec <- stepsize[2]*(1:10 - 5.5)
    temp <- matrix(NA, nrow = 4, ncol = 100)
    temp[1, ] <- rep(b1vec, each = 10)
    temp[2, ] <- rep(lamvec, 10)

    for (i in 1:10) for (k in 1:10) {
        ind <- (i - 1)*10 + k
        mtemp <- CoxEL(y=y, d=d, Z=Z, beta=b1vec[i], lam=lamvec[k], fun=Hfun)
        temp[3, ind] <- mtemp$mu
        temp[4, ind] <- mtemp$logEmpLik
    }

    subsetInd <- (temp[4, ] >= target)
    ParaTrials <- temp[, subsetInd]
    TrialValues <- Efun(beta=ParaTrials[1,], theta=ParaTrials[3,])
    maxInd <- which.max(TrialValues)
    EfunU <- TrialValues[maxInd]
    maxPara <- ParaTrials[, maxInd]
    print(c(EfunU, maxPara))

    for (N in 1:30) {
        if (sum(subsetInd) < 50) 
            stepsize <- stepsize/4
        b1vec <- maxPara[1] + stepsize[1] * (1:10 - 5.5)
        lamvec <- maxPara[2] + stepsize[2] * (1:10 - 5.5)
        temp[1, ] <- rep(b1vec, each = 10)
        temp[2, ] <- rep(lamvec, 10)
        for (i in 1:10) for (k in 1:10) {
            ind <- (i - 1) * 10 + k
            mtemp <- CoxEL(y=y, d=d, Z=Z, beta=c(b1vec[i]), lam=lamvec[k], fun=Hfun)
            temp[3, ind] <- mtemp$mu
            temp[4, ind] <- mtemp$logEmpLik
        }

        subsetInd <- (temp[4, ] >= target)
        ParaTrials <- temp[, subsetInd]
        TrialValues <- Efun(beta=ParaTrials[1,], theta=ParaTrials[3,])
        maxInd <- which.max(TrialValues)
        EfunU <- TrialValues[maxInd]
        maxPara <- ParaTrials[, maxInd]
        print(c(EfunU, maxPara))
    }
    list(Upper = EfunU, maxParameterNloglik = maxPara)
}