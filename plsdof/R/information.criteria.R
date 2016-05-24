information.criteria=function (RSS, DoF, yhat=NULL, sigmahat, n,criterion="bic"){
    if (criterion=="aic"){
        score <- as.vector(RSS/n + 2 * (DoF/n) * sigmahat^2)
    }
    if (criterion=="bic"){
        score <- as.vector(RSS/n + log(n) * (DoF/n) * sigmahat^2)
}
        
    if (criterion=="gmdl"){
       SS<-sigmahat^2
        denominator<-DoF*SS
        FF <- (yhat)/(DoF * SS)
        FF[1,FF==0]=Inf
        score <- as.vector((n/2) * log(SS) + (DoF/2) * log(FF) + (1/2) * log(n))
}
    par<-first.local.minimum(score)

    return(list(DoF = DoF, par = par, score=score))
}
