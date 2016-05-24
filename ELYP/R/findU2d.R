findU2d <- function(NPmle, ConfInt, LogLikfn, Pfun, dataMat, level=3.84)
{
#### If saving memory space is not the top priority, and search 10 by 10 trials are OK, this is the function.
#### NPmle:  should be a vector of length 2, they are the NPmle of beta1 and beta2.
#### ConfInt:  should be a vector of length 2 ... the aprox. length of two 90% confidence intervals
#### for beta1, beta2. They do not have to be exact. But they provides good initial step.
#### LogLikfn:  is a function that takes two inputs: para + dataMat, and returns the loglik value.
#### Pfun:  is a function of 2 NPmle. Pfun(beta1, beta2) = parameter we want to find confidence Interval.
#### dataMat:  is only used for the function LogLikfn( ).
####
#### Output:  The beta1 and beta2, and Pfun which is maxed inside 95 percent contour. and Loglik value.

    temp0 <- LogLikfn(mle=NPmle, dataMat)
    MaxValue <- temp0$Loglik
    PfunUold <- Pfun(NPmle[1], NPmle[2])  
    stepsize <- ConfInt/3
    b1vec <- NPmle[1] + stepsize[1] * (1:10 - 5.5)   ### pertubations, 10 by 10 or others.
    b2vec <- NPmle[2] + stepsize[2] * (1:10 - 5.5)
    
    temp <- matrix(NA, nrow = 4, ncol = 100)
    temp[1, ] <- rep(b1vec, each=10)
    temp[2, ] <- rep(b2vec, times=10)
    
    for (i in 1:10) for (j in 1:10) {
         ind <- (i - 1) * 10 + j
         temp[4, ind] <- LogLikfn(mle=c(b1vec[i],b2vec[j]),dataMat)$Loglik
         #### temp[3, ind] <- Pfun(b1vec[i], bevec[j])         #### Or other fun??
    }

    subsetInd <- (temp[4, ] >= (MaxValue - level/2))        
    ParaTrials <- temp[, subsetInd]
    TrialValues <- Pfun( ParaTrials[1,], ParaTrials[2,] ) 
    #### TrialValues <- ParaTrials[3,]

    maxInd <- which.max(TrialValues)
    PfunU <- TrialValues[maxInd]
    maxParaOLD <- maxPara <- ParaTrials[, maxInd]     
    print(c(PfunU, maxPara))

    for (N in 1:25){                #### may be more steps?
        if( sum(subsetInd) < 20) stepsize <- stepsize/2
        #### if( abs(Ub1 - mean(b1vec)) < 4*min(stepsize) ) stepsize <- stepsize/2  ####
        if(PfunU <= PfunUold) {maxpara <- maxParaOLD ; stepsize <- stepsize/4}
        b1vec <- maxPara[1] + stepsize[1]*(1:10 - 5.5)
        b2vec <- maxPara[2] + stepsize[2]*(1:10 - 5.5)
        
        temp[1, ] <- rep(b1vec, each=10)
        temp[2, ] <- rep(b2vec, times=10)
        
     for (i in 1:10) for (j in 1:10) {
            ind <- (i-1)*10 + j
            temp[4, ind] <- LogLikfn(mle=c(b1vec[i], b2vec[j]), dataMat)$Loglik
            #### temp[3, ind] <- Pfun(b1vec[i], b2vec[j])
        }
        subsetInd <- (temp[4, ] >= (MaxValue - level/2))
        ParaTrials <- temp[, subsetInd]
        TrialValues <- Pfun( ParaTrials[1,], ParaTrials[2,] )
        #### TrialValues <- ParaTrials[3,]
        maxInd <- which.max(TrialValues)
        PfunUold <- PfunU
        PfunU <- TrialValues[maxInd]
        maxPara <- ParaTrials[, maxInd]
        print(c(PfunU, maxPara))
    }
    list(Upper = PfunU, maxParameterNloglik = maxPara)
}