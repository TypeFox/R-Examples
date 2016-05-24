BJfindL2 <- function(NPmle, ConfInt, LLRfn, Betafun, dataMat, level=3.84) 
{
    if (length(NPmle) != 2) 
            stop("NPmle must of length 2")
    #### temp0 <- LLRfn(para=NPmle, dataMat)
    #### MaxValue <- temp0$"-2LLR"          #### No need for Max    
    PfunLold <- Betafun(b1=NPmle[1],b2=NPmle[2]) 
   
    stepsize <- ConfInt/4
    b1vec <- NPmle[1] + stepsize[1] * (1:6 - 3.5)
    b2vec <- NPmle[2] + stepsize[2] * (1:6 - 3.5)
    temp <- matrix(NA, nrow = 3, ncol = 36)
    temp[1, ] <- rep(b1vec, each = 6)
    temp[2, ] <- rep(b2vec, 6)

    for (i in 1:6) for (j in 1:6) {
        ind <- (i - 1) * 6 + j
        mtemp <- LLRfn(para=c(b1vec[i], b2vec[j]), dataMat)
        temp[3, ind] <- mtemp$"-2LLR"
    }
    
    target <- level
    subsetInd <- (temp[3, ] <= target ) 
    if( !any(subsetInd) ) stop("ConfInt too wide")     ##what if subsetInd is empty??
    ParaTrials <- temp[, subsetInd]
    TrialValues <- Betafun(b1=ParaTrials[1,], b2=ParaTrials[2,])
    minInd <- which.min(TrialValues)
    PfunL <- TrialValues[minInd]
    minParaOLD <- minPara <- ParaTrials[, minInd]
    
    for (N in 1:25) {
        if (sum(as.numeric(subsetInd)) < 15) 
           stepsize <- stepsize/2
        if ( PfunL >= PfunLold ) 
              { minPara <- minParaOLD ; stepsize <- stepsize/3 }
        b1vec <- minPara[1] + stepsize[1] * (1:6 - 3.5)
        b2vec <- minPara[2] + stepsize[2] * (1:6 - 3.5)
    
        temp[1, ] <- rep(b1vec, each = 6)
        temp[2, ] <- rep(b2vec, 6)
      
        for (i in 1:6) for (j in 1:6) {
            ind <- (i - 1) * 6 + j
            mtemp <- LLRfn(para = c(b1vec[i], b2vec[j]), dataMat)
            temp[3, ind] <- mtemp$"-2LLR"
        }
        subsetInd <- (temp[3, ] <= target)
        ParaTrials <- temp[, subsetInd]
        TrialValues <- Betafun(b1=ParaTrials[1,], b2=ParaTrials[2,])
        minInd <- which.min(TrialValues)
        PfunLold <- PfunL
        PfunL <- TrialValues[minInd]
        minParaOLD <- minPara
        minPara <- ParaTrials[, minInd]
        print(c(PfunL, minPara))
    }
    list(Lower = PfunL, minParameterNloglik = minPara)
}
