findL3 <- function(NPmle, ConfInt, LogLikfn, Pfun, level=3.84, dataMat)
{
#### This is the sister funtion of findU3( ), only for the lower confidence bound.
    MaxValue <- LogLikfn(c(NPmle, 0), dataMat)$Loglik
    stepsize <- ConfInt/3
    b1vec <- NPmle[1] + stepsize[1] * (1:10 - 5.5)
    b2vec <- NPmle[2] + stepsize[2] * (1:10 - 5.5)
    lamvec <- stepsize[3] * (1:10 - 5.5)
    temp <- matrix(NA, nrow = 5, ncol = 1000)
    temp[1, ] <- rep(b1vec, each=100)
    temp[2, ] <- rep(rep(b2vec, each=10), 10)
    temp[3, ] <- rep(lamvec, 100)
    for (i in 1:10) for (j in 1:10) for (k in 1:10) {
        ind <- (i-1)*100 + (j-1)*10 + k
        mtemp <- LogLikfn(mle=c(b1vec[i],b2vec[j],lamvec[k]),dataMat)
        temp[4, ind] <- mtemp$Mulam
        temp[5, ind] <- mtemp$Loglik
    }
    subsetInd <- (temp[5, ] >= (MaxValue - level/2))
    ParaTrials <- temp[, subsetInd]
    TrialValues <- Pfun(b1=ParaTrials[1,], b2=ParaTrials[2,], Mulam=ParaTrials[4,])
    minInd <- which.min(TrialValues)
    PfunL <- TrialValues[minInd]
    minPara <- ParaTrials[, minInd]
    print(c(PfunL, minPara))
    for (N in 1:18) {
        if( sum(subsetInd) < 800 ) stepsize <- stepsize/2
        b1vec <- minPara[1] + stepsize[1]*(1:10 - 5.5)
        b2vec <- minPara[2] + stepsize[2]*(1:10 - 5.5)
        lamvec <- minPara[3] + stepsize[3]*(1:10 -5.5)
        temp[1, ] <- rep(b1vec, each=100)
        temp[2, ] <- rep( rep(b2vec, each=10), 10)
        temp[3, ] <- rep(lamvec, 100)
     for (i in 1:10) for (j in 1:10) for (k in 1:10) {
            ind <- (i-1)*100 + (j-1)*10 +k
            mtemp <- LogLikfn(mle=c(b1vec[i], b2vec[j],lamvec[k]), dataMat)
            temp[4, ind] <- mtemp$Mulam
            temp[5, ind] <- mtemp$Loglik 
        }
        subsetInd <- (temp[5, ] >= (MaxValue - level/2))
        ParaTrials <- temp[, subsetInd]
        TrialValues <- Pfun(b1=ParaTrials[1,], b2=ParaTrials[2,], Mulam=ParaTrials[4,])
        minInd <- which.min(TrialValues)
        PfunL <- TrialValues[minInd]
        minPara <- ParaTrials[, minInd]
        print(c(PfunL, minPara))
    }
    list(Lower = PfunL, minParameterNloglik = minPara)
}
