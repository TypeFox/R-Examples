readPriorFromFile=function(priorFileName) {
    fPrior=file(priorFileName, "r")
    componentCnt=scan(fPrior, n=1, quiet=T)
    AA_KINDS=scan(fPrior, n=1, quiet=T)
    mixtureCoef=numeric(componentCnt)
    mAlpha=matrix(0, nrow=componentCnt, ncol=AA_KINDS);
    for (i in  1:componentCnt){
        mixtureCoef[i]=scan(fPrior, n=1, quiet=T)
        mAlpha[i,]=scan(fPrior, n=AA_KINDS, quiet=T)
    }
    prior=list(alpha=mAlpha, mix.coef=mixtureCoef)
    attr(prior, "label")="9 clouds"
    close(fPrior)
    prior
}


hmmMargLlik=function(dat,aaPrior, tau) {
    aaCounts=alignment2count(dat, level=20)
    tranCounts=alignment2trancount(dat, weight=rep(1,nrow(dat)))
    loglik=apply(aaCounts, 1, function (y) logIntegrateMixDirichlet(y, aaPrior, tau))
    loglik=sum(loglik) # summing over the 9 components in the cloud 9 prior
    loglik=loglik + sum(apply(tranCounts[,1:2,drop=FALSE], 1, function (y) logIntegrateDirichlet(y, c(.7939, .0135))))
    loglik=loglik + sum(apply(tranCounts[,3:4,drop=FALSE], 1, function (y) logIntegrateDirichlet(y, c(.9002, .5630))))
    return(loglik)
}
