### calculate a whole path for the exact results
BMNExact = function(X, rhoVec, thrCheck=1e-3, thrPseudo=1e-5, ThetaStart=NULL, verbose = FALSE, maxIter=100, timeout=60, penalize.diag=FALSE)
{
    ### initialize variables
    rhoVec = sort(rhoVec, decreasing=T)
    n = dim(X)[1]
    p = dim(X)[2]

    rhoLen = length(rhoVec)
    ThetaRes = vector("list", length=rhoLen)
    success = logical(rhoLen)
    Theta0=ThetaStart
    completed = 0
    
    for(i in 1:rhoLen)
    {
        if(verbose)
        {   print(paste("Rho=", rhoVec[i]))}
        
        res = try(BMNExact.single(X,rhoVec[i], thrCheck, thrPseudo, Theta0, verbose, maxIter, timeout, penalize.diag))
        if(class(res)=="try-error")
        {
            break
        }
        success[i]=res$success
        Theta0=res$Theta
        ThetaRes[[i]]=res$Theta
        completed = completed + 1
    }
    return(list(rho=rhoVec[1:completed],ThetaList=ThetaRes[1:completed], success=success[1:completed], penalize.diag=penalize.diag))
}

### write a program that only evaluates the whole covariance matrix of the ising model when necessary
BMNExact.single = function(X, rho, thrCheck=1e-3, thrPseudo=1e-5, ThetaStart=NULL, verbose = FALSE, maxIter=100, timeout=60, penalize.diag=FALSE)
{
    ### calculate some parameters
    n=dim(X)[1]
    p=dim(X)[2]
    S=t(X) %*% X/n

    ### calculate some needed values and process input parameters
    if(is.null(ThetaStart))
    {
        ThetaStart = matrix(numeric(p^2), ncol=p)
        diag(ThetaStart) = log(diag(S)/(1-diag(S)))
        ### dont let it be too large
        ThetaStart[abs(ThetaStart) > 1] = 1
    }
    Theta=ThetaStart
    rhoMat = BMNcheckRho(rho,p,penalize.diag)
    converged=FALSE
    success=FALSE
    curTime=proc.time()[1]
    
    runAll = try({
    W = BMNJT(Theta, adjMat=NULL, var=NULL, onlyActive=FALSE, (proc.time()[1]-curTime))$SecondMomentMatrix
    iterOuter = 1
    
    ### start the outer loop; always uses the full covariance matrix
    while(iterOuter <= maxIter && !converged)
    {
        ### use the full covariance matrix to make one step with the pseudo function, except in the first step
        if(iterOuter==1){Delta = matrix(numeric(p^2), ncol=p)}else{Delta = BMNcalcDelta(X, Theta, W)}
        res = BMNPseudo.single(X,rhoMat,Delta, ThetaStart=Theta, maxError=thrPseudo, penalize.diag=penalize.diag)
        if(!res$success)
        {
            break
        }
        else
        {
            Theta = BMNgetNextTheta(res$Theta, Theta, S, W, rhoMat)
        }
        ### inner loop only uses non-zero variables in the optimization. All others are forced to be 0
        iterInner = 1
        convergedInner = FALSE
        ### save non-zero elements to always optimize over the same variables in the inner loop
        adjMat = (Theta!=0)
        rhoMatRest = BMNgetRestricedRhoMatrix(Theta, rhoMat)
        W = BMNJT(Theta, adjMat=adjMat, var=NULL, onlyActive=TRUE, (proc.time()[1]-curTime))$SecondMomentMatrix
        while(iterInner <= maxIter && !convergedInner)
        {
            Delta = BMNcalcDelta(X,Theta,W)
            res = BMNPseudo.single(X,rhoMatRest,Delta, ThetaStart=Theta, maxError=thrPseudo, penalize.diag=penalize.diag)
            if(!res$success)
            {
                stop("Pseudo-likelihood failed")
            }
            else
            {
                Theta = BMNgetNextTheta(res$Theta, Theta, S, W, rhoMatRest)
            }
            W = BMNJT(Theta, adjMat=adjMat, var=NULL, onlyActive=TRUE, (proc.time()[1]-curTime))$SecondMomentMatrix
            ### check for convergence
            res=BMNcheckSubGrad(W,S,Theta,rhoMatRest, thrCheck)
            if(verbose)
            {
                print(paste("Inner Iteration:",iterInner))
                print(res)
            }
            if(res$converged)
            {
                convergedInner=TRUE
            }
            else
            {
                Delta = BMNcalcDelta(X, Theta, W)
            }
            if(proc.time()[1]-curTime>timeout)
            {
                stop("Timeout")
            }
            iterInner = iterInner + 1
        }
        ### check for convergence in the outer loop
        if(iterInner > maxIter) ### inner loop did not converge; stop with failure
        {
            stop("Too many iterations in inner loop")
        }
        W = BMNJT(Theta, adjMat=NULL, var=NULL, onlyActive=FALSE, (proc.time()[1]-curTime))$SecondMomentMatrix
        ### check for convergence
        res=BMNcheckSubGrad(W,S,Theta,rhoMat, thrCheck)
        if(verbose)
        {
            print(paste("Outer Iteration:",iterOuter))
            print(res)
        }
        if(res$converged)
        {
            converged=TRUE
            success=TRUE
        }
        else
        {
            Delta = BMNcalcDelta(X, Theta, W)
        }
        
        if(proc.time()[1]-curTime>timeout)
        {
            stop("Timeout")
        }

        iterOuter = iterOuter+1
    }
    }, silent=T)
    if(class(runAll)=="try-error")
    {
        success=FALSE
        Theta=NULL
    }
    
    return(list(Theta=Theta, success=success))
}



BMNgetNextTheta = function(newTheta, oldTheta, S, W, rhoMat)
{
    ### do a line search step
    gradient = BMNGradientPen(W,S,oldTheta, rhoMat)
    dTheta = newTheta-oldTheta
    stepSize = BMNlineSearch(oldTheta, S, gradient, dTheta, rhoMat)
    newTheta=oldTheta + stepSize*dTheta
    ### throw away small values
    newTheta = newTheta * (abs(newTheta) > 1e-6)
    return(newTheta)
}

### creates a rho matrix that has penalty infinity wherever Theta is 0
BMNgetRestricedRhoMatrix = function(Theta, rhoMat)
{
    rhoMatRest = rhoMat
    rhoMatRest[Theta==0] = 10^10 ## considered to be infinity
    return(rhoMatRest)
}


BMNcheckRho = function(rho,p, penalize.diag)
{
    if(is.matrix(rho)) # check that it is a matrix
    {
        if(isTRUE(dim(rho)==c(p,p))) # square
        {
            if(isTRUE(rho>=0)) # with non-negative elements
            {
                return(rho)
            }
            else
            {
                stop("entries of rho have to be non-negative")
            }
        }
        else
        {
            stop("rho has to be a p by p matrix or a scalar")
        }
    }
    else # check that it is a scalar
    {
        if(length(rho)==1)
        {
            if(rho[1]>=0)
            {
                rhoMat = matrix(numeric(p^2), ncol=p)
                rhoMat[,]=rho
                if(!penalize.diag)
                {
                    diag(rhoMat)=0
                }
                return(rhoMat)
            }
            else
            {
                stop("rho has to be non-negative")
            }
        }
        else
        {
            stop("rho has to be a p by p matrix or a scalar")
        }
    }
}

BMNcalcDelta = function(X, Theta, W)
{
    n = dim(X)[1]
    ### calculate the pseudo-likelihood probabilities and the adjustment Delta
    etaMat=BMNcalcEtas(X,Theta)
    pMat = 1-1/(1+exp(etaMat))
    foo = t(pMat) %*% X
    diag(foo) = apply(pMat, MARGIN=2, FUN=sum)
    foo = foo + t(foo)
    diag(foo) = apply(pMat, MARGIN=2, FUN=sum)
    diag(foo) = diag(foo) + apply(X^2,MARGIN=2,FUN=sum)
    bar = 2*n*W
    diag(bar) = diag(bar)
    Delta = bar-foo
    return(Delta)
}


BMNcalcEtas = function(X,Theta)
{
    n = dim(X)[1]
    p = dim(X)[2]
    etaMat = matrix(numeric(n*p), ncol=p)
    for(j in 1:p)
    {
        etaMat[,j] = BMNcalcSingleEta(X,Theta,j)
    }
    return(etaMat)
}

### calculate a single probability
BMNcalcSingleEta = function(X,Theta,j)
{
    Xtilde = X
    Xtilde[,j]=1
    return(Xtilde %*% Theta[,j])
}



BMNlineSearch=function(Theta, S, gradient, descDirec, rhoMat, timeout, alpha=0.1,beta=0.8,maxIter=20)
{
    # check that it is a descent direction
    if(sum(gradient*descDirec)>0) # not a descent direction
    {   t= -1}
    else
    {   t=1}
    
    # do backtracking
    newTheta=Theta+t*descDirec
    fOld = BMNL1Loss(S,Theta,rhoMat)
    fNew = BMNL1Loss(S,newTheta,rhoMat)
    iter=1
    while(fNew > fOld + alpha*t*sum(gradient*descDirec) && iter<maxIter)
    {
        iter=iter+1
        t=t*beta
        newTheta=Theta+t*descDirec
        fNew = BMNL1Loss(S,newTheta,rhoMat)
    }
    
    return(t)
}


### function that copies the upper right triangle of a matrix onto the lower left
BMNmakeSymmetric=function(mat)
{
    mat[lower.tri(mat)] = (t(mat))[lower.tri(mat)]
    return(mat)
}


### calculate the value of the loss function
BMNL1Loss = function(S, Theta, rhoMat)
{
    Theta=BMNmakeSymmetric(Theta)
    logPartFunc = BMNJTlogPartFunc(Theta)
    Theta[lower.tri(Theta)]=0
    return(logPartFunc-sum(S*Theta)+sum(abs(rhoMat*Theta)))
}

### calculate the value of the logLikelihood
BMNLogLik = function(S, Theta)
{
    Theta=BMNmakeSymmetric(Theta)
    logPartFunc = BMNJTlogPartFunc(Theta)
    Theta[lower.tri(Theta)]=0
    return(-logPartFunc+sum(S*Theta))
}


### calculate the value of gradient paying special attention to the case when
### an active variable is 0
BMNGradientPen = function(W,S,Theta, rhoMat)
{
    gradient = (W-S+rhoMat*sign(Theta))
    Zero = (Theta==0)
    gradZero = (W-S)[Zero]
    rhoMatZero = rhoMat[Zero]
    gradZero[abs(gradZero)<rhoMat[Zero]]=0
    foo = abs(gradZero)>=rhoMat[Zero]
    gradZero[foo]=gradZero[foo]-sign(gradZero[foo])*rhoMatZero[foo]
    gradient[Zero]=gradZero
    return(gradient)
}



### function that checks if the subgradient equations are satisfied
BMNcheckSubGrad = function(W,S, Theta, rhoMat, err=0.001)
{
    Zero = (rhoMat == 0)
    rhoMat[rhoMat== 0] = err
    T = (S-W)/rhoMat
    violLarge = sum(abs(T)>1+err)
    violSign = sum(((sign(Theta) != sign(T)) & abs(Theta)>err) & !Zero)
    violT1 =  sum(abs(Theta)>err & abs(T)<1-err & !Zero) 
    browser()
    if((violLarge==0) && (violSign==0) && (violT1==0))
    {
        converged=TRUE
    }
    else
    {
        converged=FALSE
    }
    return(list(converged=converged, violLarge=violLarge, violSign=violSign, violT1=violT1))
}


