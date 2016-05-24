
MCMCfixedEF = function(nn,PP,KK,dd,XX,YY,ZZ,TT,
		beta,intercept,alpha,
		MuAlpha,SigmaAlpha,MuBeta,SigmaBeta,MuZ,VarZ,
		tuneBetaAll,tuneInt,tuneAlpha,tuneZAll,
		niter,PriorA,PriorB,intervention )
{
    accBetaAll = rep(0,PP) 
    accAlpha = 0
    accIntAll = 0
    accZAll = rep(0,sum(nn))
    betaFinal = rep(0,niter*length(beta))
    AlphaFinal = rep(0, niter*length(alpha))
    ZZFinal = rep(0,niter*length(ZZ))
    InterceptFinal = rep(0, niter*length(intercept))
    Zvar1 = Zvar2= rep(0, niter)
    likelihood = rep(0,niter)

    out = .C('sampleFixedIntervention',as.integer(niter), as.double(XX), as.double(YY), as.double(ZZ), as.integer(TT), as.integer(nn), as.integer(PP), as.integer(dd), as.integer(KK), as.double(beta), as.double(intercept), as.double(alpha), as.double(MuAlpha), as.double(SigmaAlpha), as.double(MuBeta), as.double(SigmaBeta), as.double(MuZ), as.double(VarZ), as.double(tuneBetaAll), as.double(tuneInt),as.double(tuneAlpha), as.double(tuneZAll), as.double(accBetaAll), as.double(accAlpha), as.double(accIntAll),as.double(accZAll), as.double(betaFinal), as.double(AlphaFinal), as.double(ZZFinal), as.double(InterceptFinal),as.double(Zvar1), as.double(Zvar2),as.double(likelihood),as.double(PriorA),as.double(PriorB), as.integer(intervention))

    betaFinal = array(out[[27]],dim = c(niter,PP,1))
    AlphaFinal = array(out[[28]],dim = c(niter,1) )
    InterceptFinal = array(out[[30]],dim=c(niter,1) )

    ZZFinal = list()
    accZ = list()
    ZZx = array(out[[29]],dim = c((sum(nn)*dd),niter))
    for(ni in 1:niter ){
	Zsm = list()
        n00 = n0 = 1	
    for(kk in 1:KK ){
        n1 = sum(nn[1:kk] )
        Zsm[[kk]] = array(ZZx[n0:(2*n1),ni],dim = c(nn[kk],dd))
        accZ[[kk]] = out[[26]][n00:(n1)]/(niter)
        n00 = (n1) + 1
	n0 = (2*n1) + 1
}
	ZZFinal[[ni]] = Zsm
}

#    n00 = n0 = 1
#    ZZx = array(out[[29]],dim = c(sum(nn),dd,niter))
#    for(kk in 1:KK ){
#        n1 = sum(nn[1:kk] )
#        ZZFinal[[kk]] = array(ZZx[n00:n1,,],dim = c(nn[kk],dd,niter))
#        accZ[[kk]] = out[[26]][n00:(n1)]/(niter)
#        n00 = (n1) + 1
#}
    
    Zvar1 = out[[31]]
    Zvar2 = out[[32]]
    Zvar = data.frame(Zvar1,Zvar2)
    likelihood = out[[33]]

    accbeta = out[[23]]/(niter)
    accalpha = out[[24]]/niter
    accint = out[[25]]/(niter)
    
    draws = list(Intercept = InterceptFinal, Beta= betaFinal, Alpha = AlphaFinal, ZZ = ZZFinal, Zvar = Zvar,likelihood = likelihood)

    accrate = list(intercept = accint, beta = accbeta, alpha = accalpha,Z=accZ)

    return(list(draws = draws, acc = accrate)) 
    rm(out)
}
