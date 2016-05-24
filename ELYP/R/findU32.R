findU32 <- function(NPmle, ConfInt, LogLikfn, Pfun, dataMat, level=3.84)
{
#### NPmle: should be a vector of length 2, they are the NPmle 
####      of beta1 and beta2.
#### ConfInt: should be a vector of length 3...the length of   
####        three 90% confidence intervals,
#### beta1, beta2 and lambda. It do not have to be exact. But  
####       they provides good initial step.
#### LogLikfn: is a function that takes two inputs: para and
####           dataMat, and returns the loglik value
#### Pfun: is a function that takes the input of 3 parameter
####        values (beta1,beta2 and Mulam) and 
####        returns a parameter that we wish to find conf 
####         Interval Upper Value.
#### dataMat:  is only used for the function LogLikfn( ).
####
#### Output:  Maximum value of Pfun; and corresponding values 
####     of beta1, beta2, lambda and Mulam. and Loglik value.
    MaxValue <- LogLikfn(c(NPmle, 0), dataMat)$Loglik    ## since the NPMLE of lam is 0 always
    stepsize <- ConfInt/2
    b1vec <- NPmle[1] + stepsize[1] * (1:8 - 4.5)   ### pertubations
    b2vec <- NPmle[2] + stepsize[2] * (1:8 - 4.5)
    lamvec <- stepsize[3] * (1:8 - 4.5)
    temp <- matrix(NA, nrow = 5, ncol = 512)
    temp[1, ] <- rep(b1vec, each=64)
    temp[2, ] <- rep(rep(b2vec, each=8), 8)
    temp[3, ] <- rep(lamvec, 64)
    for (i in 1:8) for (j in 1:8) for (k in 1:8) {
        ind <- ((i - 1) * 64) + ((j - 1) * 8) + k
        mtemp <- LogLikfn(mle=c(b1vec[i],b2vec[j],lamvec[k]),dataMat)
        temp[4, ind] <- mtemp$Mulam
        temp[5, ind] <- mtemp$Loglik
    }
    subsetInd <- (temp[5, ] >= (MaxValue - level/2))
    ParaTrials <- temp[, subsetInd]
    TrialValues <- Pfun(b1=ParaTrials[1,], b2=ParaTrials[2,], Mulam=ParaTrials[4,])
    maxInd <- which.max(TrialValues)
    PfunU <- TrialValues[maxInd]
    maxPara <- ParaTrials[, maxInd]
    print(c(PfunU, maxPara))
    for (N in 1:25) {                ## may be more steps?
        if( sum(subsetInd) < 330 ) stepsize <- stepsize/1.6        ## may be do not shorten the step as fast?
        b1vec <- maxPara[1] + stepsize[1]*(1:8 - 4.5)
        b2vec <- maxPara[2] + stepsize[2]*(1:8 - 4.5)
        lamvec <- maxPara[3] + stepsize[3]*(1:8 - 4.5)
        temp[1, ] <- rep(b1vec, each=64)
        temp[2, ] <- rep( rep(b2vec, each=8), 8)
        temp[3, ] <- rep(lamvec, 64)
     for (i in 1:8) for (j in 1:8) for (k in 1:8) {
            ind <- (i-1)*64 + (j-1)*8 +k
            mtemp <- LogLikfn(mle=c(b1vec[i], b2vec[j],lamvec[k]), dataMat)
            temp[4, ind] <- mtemp$Mulam
            temp[5, ind] <- mtemp$Loglik 
        }
        subsetInd <- (temp[5, ] >= (MaxValue - level/2))
        ParaTrials <- temp[, subsetInd]
        TrialValues <- Pfun(b1=ParaTrials[1,], b2=ParaTrials[2,], Mulam=ParaTrials[4,])
        maxInd <- which.max(TrialValues)
        PfunU <- TrialValues[maxInd]
        maxPara <- ParaTrials[, maxInd]
        print(c(PfunU, maxPara))
    }
    list(Upper = PfunU, maxParameterNloglik = maxPara)
}