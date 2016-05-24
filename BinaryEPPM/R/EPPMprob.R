EPPMprob <-
function (vlambda) {
    nmax1 <- length(vlambda)
    nmax  <- nmax1 - 1
    vwork <- vlambda[1:nmax]
    tran  <- cbind(rbind(diag(-vlambda, nmax1))) +
             cbind(0,rbind(diag(vwork, nmax),0))
    tran  <- expm(tran)
    prob  <- tran[1,]
    return(prob)                }
