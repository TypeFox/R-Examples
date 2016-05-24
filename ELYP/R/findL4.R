findL4 <- function(NPmle, ConfInt, LogLikfn2, Pfun, dataMat, level=3.84)
{
#### NPmle:  should be a vector of length 3, they are the NPmle of beta1, beta2, and a=alpha.
#### ConfInt:  should be a vector of length 4.... the length of four 90% confidence intervals, for
#### beta1, beta2, a and lambda. It do not have to be exact. But they provides good initial step.
#### LogLikfn2:  is a function that takes two inputs: para + dataMat, and returns the loglik value
#### Pfun:  is a function that takes the input of 3 parameter values (beta1,beta2,and Mulam) and 
####         returns a parameter that we wish to find conf Interval Upper Value. 
#### dataMat:  is only used for the function LogLikfn2( ).
####
#### Output:  Max Pfun; corresponding beta1, beta2, lambda and Mulam. and Loglik value.
##########################################################################################
#### Some possible improvement: shorten the stepsize diff in diff direction. How to determine
####        if stepsize is too large on one direction?

    if(length(NPmle) != 3) stop("NPmle must of length 3")
    temp0 <- LogLikfn2(c(NPmle, 0), dataMat)    ## since the NPMLE of lam is 0 always
    MaxValue <- temp0$Loglik
    PfunLold <- Pfun(b1=NPmle[1], b2=NPmle[2], Mulam=temp0$Mulam)
    print(c(PfunLold, MaxValue - level/2))
    stepsize <- ConfInt/2
    b1vec <- NPmle[1] + stepsize[1] * (1:6 - 3.5)   ### pertubations
    b2vec <- NPmle[2] + stepsize[2] * (1:6 - 3.5)
    avec <- NPmle[3] + stepsize[3] * (1:6 - 3.5)
    lamvec <- stepsize[4] * (1:6 - 3.5)

    temp <- matrix(NA, nrow = 6, ncol = 1296)   ## 8*8*8*8 = 512*8 = 4096
    temp[1, ] <- rep(b1vec, each=216)           ## 8*8 = 64, 8*8*8=512
    temp[2, ] <- rep(rep(b2vec, each=36), 6)
    temp[3, ] <- rep(rep(avec, each=6), 36)
    temp[4, ] <- rep(lamvec, 216)
    for (i in 1:6) for (j in 1:6) for (k in 1:6) for (r in 1:6) {
        ind <- (i - 1)*216 + (j - 1)*36 + (k-1)*6 + r
        mtemp <- LogLikfn2(mle=c(b1vec[i],b2vec[j],avec[k],lamvec[r]),dataMat)
        temp[5, ind] <- mtemp$Mulam
        temp[6, ind] <- mtemp$Loglik
    }

    subsetInd <- (temp[6, ] >= (MaxValue - level/2))
    ParaTrials <- temp[, subsetInd]
    TrialValues <- Pfun(b1=ParaTrials[1,], b2=ParaTrials[2,], Mulam=ParaTrials[5,])
    minInd <- which.min(TrialValues)
    PfunL <- TrialValues[minInd]
    #### Pfunstep <- 3
    minParaOLD <- minPara <- ParaTrials[, minInd]
    print(c(PfunL, minPara))

    for (N in 1:20) {                #### may be more steps?
        if( sum(subsetInd) < 300 ) stepsize <- stepsize/1.6        #### may be do not shorten the step as fast?
        ####if( abs(PfunU - PfunUold) < Pfunstep )  stepsize <- stepsize/2
        ####Pfunstep <- abs(PfunUold - PfunU)
        if(PfunL >= PfunLold) {minPara <- minParaOLD ; stepsize <- stepsize/3}
        b1vec <- minPara[1] + stepsize[1]*(1:6 - 3.5)
        b2vec <- minPara[2] + stepsize[2]*(1:6 - 3.5)
        avec <- minPara[3] + stepsize[3]*(1:6 - 3.5)
        lamvec <- minPara[4] + stepsize[4]*(1:6 - 3.5)
        temp[1, ] <- rep(b1vec, each=216)
        temp[2, ] <- rep(rep(b2vec, each=36), 6)
        temp[3, ] <- rep(rep(avec, each=6), 36)
        temp[4, ] <- rep(lamvec, 216)
     for (i in 1:6) for (j in 1:6) for (k in 1:6) for (r in 1:6){
            ind <- (i-1)*216 + (j-1)*36 + (k-1)*6 + r
            mtemp <- LogLikfn2(mle=c(b1vec[i], b2vec[j], avec[k], lamvec[r]), dataMat)
            temp[5, ind] <- mtemp$Mulam
            temp[6, ind] <- mtemp$Loglik 
        }
        subsetInd <- (temp[6, ] >= (MaxValue - level/2))
        ParaTrials <- temp[, subsetInd]
        TrialValues <- Pfun(b1=ParaTrials[1,], b2=ParaTrials[2,], Mulam=ParaTrials[5,])
        minInd <- which.min(TrialValues)
        PfunLold <- PfunL
        PfunL <- TrialValues[minInd]
        minParaOLD <- minPara
        minPara <- ParaTrials[, minInd]
        print(c(PfunL, minPara))
    }
    list(Lower = PfunL, minParameterNloglik = minPara)
}
