BMNPseudo.single = function(X, rho, Delta=NULL, ThetaStart=NULL, maxError=1e-5, maxIter=100, penalize.diag=FALSE, stepSize=1, performLineSearch=FALSE)
{
    ### parameters
    n= dim(X)[1]
    p = dim(X)[2]
    XTX = t(X) %*% X
    
    ### check that there is at least one observation 0 or 1 in each variable
    foo = diag(XTX)
    if(min(foo)==0 || max(foo)==n)
    {
        stop("Each variable has to have at least one observation that is 0 and 1.")
    }
    
    ### set rho right only if it is not a matrix
    if(is.matrix(rho))
    {
        R_rho = 2*n*rho
    }
    else
    {
        R_rho = matrix(numeric(p^2), ncol=p)
        R_rho[,] = 2*n*rho
        if(!penalize.diag){diag(R_rho)=0}
    }
    
    ### set Delta
    if(is.null(Delta)){Delta = matrix(numeric(p^2), ncol=p)}
    ### set ThetaStart
    if(is.null(ThetaStart))
    {
        ThetaStart = matrix(numeric(p^2), ncol=p)
        diag(ThetaStart) = log(diag(XTX)/(n-diag(XTX)))
        ### dont let it be too large
        ThetaStart[abs(ThetaStart) > 1] = 1
    }
    
    ### iteration numbers and other settings
    maxOuterIter = as.integer(max(maxIter,p))
    maxInnerIter = as.integer(10000)
    maxVarAdd = as.integer(4*p)
    performLineSearch = as.logical(performLineSearch)
    
    ### run the analysis
    foo=.Call("pseudoLikelihood",X, XTX, R_rho, Delta, ThetaStart, maxError, maxOuterIter, maxInnerIter, maxVarAdd, stepSize, performLineSearch)

    return(foo)
}


BMNPseudo = function(X, rhoVec, Delta=NULL, ThetaStart=NULL, maxError=1e-5, verbose=FALSE, maxIter=100, penalize.diag=FALSE, stepSize=1, performLineSearch=FALSE)
{
    ### initialize variables
    rhoVec = sort(rhoVec, decreasing=T)
    n = dim(X)[1]
    p = dim(X)[2]

    rhoLen = length(rhoVec)
    ThetaRes = vector("list", length=rhoLen)
    success = logical(rhoLen)
    Theta0=ThetaStart
    
    for(i in 1:rhoLen)
    {
        if(verbose)
        {   print(paste("Rho=", rhoVec[i]))}
        res = BMNPseudo.single(X,rhoVec[i], Delta, Theta0, maxError, maxIter, penalize.diag, stepSize, performLineSearch)
        success[i]=res$success
        Theta0=res$Theta
        ThetaRes[[i]]=res$Theta
    }
    return(list(rho=rhoVec[1:length(ThetaRes)],ThetaList=ThetaRes, success=success, penalize.diag=penalize.diag))
}
