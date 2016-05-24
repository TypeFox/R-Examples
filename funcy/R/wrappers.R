#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#


fitfclustWrapper <- function(data, k, reg, regTime, funcyCtrlMbc,
                             fpcCtrl, p=2, pert=0.01){
    ##input
    dimBase <- funcyCtrlMbc@dimBase
    baseType <- funcyCtrlMbc@baseType
    epsilon <- funcyCtrlMbc@eps
    maxiter <- funcyCtrlMbc@maxit
    hard <- funcyCtrlMbc@hard
    seed <- funcyCtrlMbc@seed
    init <- funcyCtrlMbc@init
    nrep <- funcyCtrlMbc@nrep
    redDim <- funcyCtrlMbc@redDim
    
    if(reg==1){
        predFct <- fitfclust.pred
        curvepredFct <- fitfclust.curvepred
        discrimFct <- fitfclust.discrim
        plotCurvesFct <- fitfclust.plotcurves
        
    }else{
        predFct <- fitfclust.predIrreg
        curvepredFct <- fitfclust.curvepredIrreg
        discrimFct <- fitfclust.discrim
        plotCurvesFct <- fitfclust.plotcurvesIrreg
    }

    ##evaluation
    ptm <- proc.time()
    res <- fitfclust(data=data, dimBase=dimBase, h=redDim, p=p, k=k,
                     regTime=regTime,
                     epsilon=epsilon, maxiter=maxiter, pert=pert,
                     hard=hard, seed=seed, init=init, nrep=nrep, reg=reg,
                     fpcCtrl=fpcCtrl, baseType=baseType)
    sysTime <- proc.time()-ptm
    

    prms <- res$parameters
    ##class prediction
    pred <- predFct(res)
    ##curve prediction
    curvePred <- curvepredFct(res)

    ##funcyOut
    out <- new("funcyOutMbc-fitfclust")
    out@methodName <- "fitfclust"
    out@kOut <- k
    out@time <- res$grid
    out@dimBaseOut <- dimBase
    out@cluster<-pred$class.pred
    out@centers <- curvePred$meancurves
    out@props <- round(prms[[4]],4)
    out@dist2centers<-dist2centers(data, out@centers)
    out@cldist=makeClMat(out@dist2centers)
    out@calcTime <- sysTime 
    out@plotParams <- res$plotParams
    ##funcyOutMbc 
    out@groupDimBase <- rep(dim(prms$alpha)[2], k)
    out@probs<-res$vars$piigivej
    out@prms <- prms[-4]
    out@AIC <- res$aic
    out@BIC <- res$bic
    out@logLik <- res$loglik
    out@nrIter <- res$nrIter
    ##funcyOutMbc-fitfclust output
    out@fit <- res #for plotting

    return(out)
}


distclustWrapper <- function(data, k, reg, regTime, funcyCtrl,
                            fpcCtrl, method="pam")
    {
        ##evaluation
        ptm <- proc.time()
        res <- distclust(data, k, reg, regTime, funcyCtrl,
                        fpcCtrl, method=method)
        sysTime <- proc.time()-ptm

        ##funcyOut
        out <- new("funcyOut")
        out@methodName <- "distclust"
        out@kOut <- k
        out@time <- res$grid
        out@dimBaseOut= res$dimBase
        out@cluster=res$cluster
        out@centers=res$centers
        out@props=round(res$props,4)
        out@dist2centers=res$dist2centers
        out@cldist=makeClMat(out@dist2centers)
        out@calcTime=sysTime 
        out@plotParams=res$plotParams
        
        return(out)
    }


iterSubspaceWrapper <- function(data, k, reg, regTime, funcyCtrlMbc,
                            fpcCtrl, simplif=TRUE)
    {
        ##evaluation
        ptm <- proc.time()
        res <- iterSubspace(data=data, k=k, regTime=regTime,
                        reg=reg, funcyCtrlMbc=funcyCtrlMbc,
                        fpcCtrl=fpcCtrl, simplif=simplif)
        sysTime <- proc.time()-ptm
        
        ##funcyOut
        out <- new("funcyOut-iterSubspace")
        out@methodName <- "iterSubspace"
        out@kOut <- k
        out@time <- res$grid
        out@dimBaseOut <- funcyCtrlMbc@dimBase
        out@cluster<-res$cls
        out@centers<-res$ctrs
        out@props <- round(as.numeric(table(res$cls)/length(res$cls)),4)
        out@dist2centers <- dist2centers(data, out@centers)
        out@cldist=makeClMat(out@dist2centers)
        out@calcTime<- sysTime 
        out@plotParams <- res$plotParams
        ##funcyOut-iterSubspace 
        out@groupDimBase<-res$groupDimBase
        out@prms <- list(groupMeans=res$groupMeans,
                         groupBase=res$groupBase,
                         groupErr=res$groupErr)
        out@nrIter <- res$nrIter

        return(out)
}


funclustWrapper <- function(data, k, reg, regTime, funcyCtrlMbc,
                            nbInit=5, nbIterInit=50, ...){

    if(!requireNamespace("Funclustering"))
        stop("Please install Funclustering to use this method.")
    
    if(!reg)
        stop("This method does not work on sparse data!")
    baseType <- funcyCtrlMbc@baseType
    if(baseType=="eigenbasis")
        stop("This base type is not implemented yet.")
    if(!is.null(funcyCtrlMbc@seed))
        warning("It is not possible to set a seed for method funclust.")

    ##input
    dimBase <- funcyCtrlMbc@dimBase
    thd <- funcyCtrlMbc@thd
    epsilon <- funcyCtrlMbc@eps
    nbIteration <- funcyCtrlMbc@maxit
    hard <- funcyCtrlMbc@hard
    fixedDimension <- rep(funcyCtrlMbc@redDim, k)
    increaseDimension <- funcyCtrlMbc@flexDim

    ##evaluation
    ptm <- proc.time()
    res <- formatFuncy(data, regTime=regTime,  format="Format3")
    data <- t(res$Yin); t_all <- res$t_all
    baseObj <- makeBasis(baseType, t_all, nbasis=dimBase)$bObj
    fd <- Data2fd(data, argvals=t_all, basisobj=baseObj);
    res=Funclustering::funclust(fd=fd, K=k, thd=thd, increaseDimension=increaseDimension, hard=hard, fixedDimension=fixedDimension, nbInit=nbInit,
        nbIterInit=nbIterInit, nbIteration=nbIteration,
        epsilon=epsilon, ...)
    sysTime <- proc.time()-ptm

    clout <- label2lowerk(res$cls)
    
    ##funcyOut
    out <- new("funcyOutMbc")
    out@methodName <- "funclust"
    out@kOut <- clout$k
    out@time <- t_all
    out@dimBaseOut <- dimBase
    out@cluster <- clout$cluster
    centers <- sapply(res$meanList[[1]], function(x)
        eval.fd(t_all, x))[,as.numeric(names(table(res$cls)))]
    out@centers <- cbind(centers)
    out@props <- round(as.numeric(res$proportions),4)
    out@dist2centers <- dist2centers(data, out@centers)
    out@cldist=makeClMat(out@dist2centers)
    out@calcTime <- sysTime 
    ##funcyOutMbc
    out@groupDimBase <- res$dimensions
    out@probs<-res$tik
    out@prms <- list(NA)
    out@AIC <- res$aic
    out@BIC <- res$bic
    out@logLik <- -res$loglikelihood
    out@nrIter <- as.integer(NA)
    
    return(out)
}


funHDDCWrapper <- function(data, k, reg, regTime,  funcyCtrlMbc,
                           model="AkBkQkDk", ...){

     if(!requireNamespace("funHDDC"))
        stop("Please install funHDDC to use this method.")
     if(!reg)
        stop("This method does not work on sparse data!")
    baseType=funcyCtrlMbc@baseType
    if(baseType=="eigenbasis")
        stop("This base type is not implemented yet.")
    if(k==0)
        stop("Clustering for this data not possible.")

    ##input
    dimBase <- funcyCtrlMbc@dimBase
    thd <- funcyCtrlMbc@thd
    eps <-funcyCtrlMbc@eps
    maxit <-funcyCtrlMbc@maxit
    init <-funcyCtrlMbc@init
    seed <- funcyCtrlMbc@seed

    ptm <- proc.time()
    res <- formatFuncy(data, regTime=regTime, format="Format3")
    data <- t(res$Yin); t_all <- res$t_all
    baseObj <- makeBasis(baseType, t_all, nbasis=dimBase)$bObj
    fd <- Data2fd(data, argvals=t_all, basisobj=baseObj);
    set.seed(seed)
    res=try(funHDDC::funHDDC(fd=fd, K=k, init=init, model=model, thd=thd,
        maxit=maxit ,eps=eps, ...), silent=TRUE)
    if(class(res)=="try-error"){
        warning(paste("Clustering with", k,
                      "classes is not possible.", k-1,
                      "clusters are used!"))
        k <- k-1
        return(funHDDCWrapper(data, k, reg, regTime, funcyCtrlMbc,
                              model="AkBkQkDk", ...))
    }
    sysTime <- proc.time()-ptm

    ##funcyOut
    out <- new("funcyOutMbc")
    out@methodName <- "funHDDC"
    out@kOut <- k
    out@time <- t_all
    out@dimBaseOut <- max(res$cls)
    out@cluster <- res$cls
    fdmeans <- fd; fdmeans$coefs <- t(res$prms$m)
    out@centers <- eval.fd(t_all, fdmeans)
    out@props <- round(as.numeric(unlist(res$prms$prop)),4)
    out@dist2centers <- dist2centers(data, out@centers)
    out@cldist=makeClMat(out@dist2centers)
    out@calcTime<- sysTime 
    ##funcyOutMbc 
    out@groupDimBase <- as.numeric(res$prms$d)
    out@probs<-res$P
    out@prms <- res$prms
    out@AIC <- -res$aic
    out@BIC <- -2*res$bic
    out@logLik <- res$loglik[length(res$loglik)]
    out@nrIter <- as.integer(NA)
    
    return(out)
}



fscmWrapper <- function(data, k, reg, regTime, funcyCtrlMbc,
                        fpcCtrl, location=NULL, scale=FALSE,
                        knn=5, useCode="C", verbose=FALSE){
    if(!reg)
        stop("This method does not work on sparse data!")

    ##evaluation
    ptm <- proc.time()
    res <- fscm(data, k, reg=reg, regTime=regTime, funcyCtrlMbc,
            location=location, scale=scale, knn=knn, useCode=useCode, verbose=verbose)
    sysTime <- proc.time()-ptm

    ##funcyOut
    out <- new("funcyOutMbc-fscm")
    out@methodName <- "fscm"
    out@kOut <- res$k
    out@time <- res$grid
    out@dimBaseOut <- res$dimBase
    out@cluster <- res$cluster
    out@centers <- res$centers
    out@props <- round(as.numeric(table(out@cluster)/length(out@cluster)),4)
    out@dist2centers <- dist2centers(data, out@centers)
    out@cldist=makeClMat(out@dist2centers)
    out@calcTime <- sysTime 
    ##funcyOutMbc
    out@groupDimBase <- rep(res$dimBase,4)
    #out@probs <- NA
    out@prms <- res$prms
    out@AIC <- -res$AIC
    out@BIC <- -res$BIC
    #out@logLik <- NA
    out@nrIter <- as.integer(NA)
   
    ##specific out object
    out@location <- res$location
    out@trends <- res$trends

    return(out)

}



waveclustWrapper <- function(data, k, reg, regTime, funcyCtrlMbc,
                             gamma="group", init="SEM", plotLoglik=FALSE){
    
                     
    if(!reg)
        stop("This method does not work on sparse data!")
    baseType=funcyCtrlMbc@baseType
    if(baseType!="wavelets")
        warning("This method works only for a wavelet basis. It will be used here.")
    
    if(is.null(regTime))
        regTime <- 1:dim(data)[1]
    ##evaluation
    Y1 <- apply(t(data),1, list)
    Y <- lapply(Y1,unlist)
    CCD <- new("CClustData", Y=Y, filter.number=1)
    CCDred <- getUnionCoef(CCD)
    CCO <- new("CClustO")
   
    CCO@nbclust <- k
    CCO@Gamma2.structure <- gamma
    CCO@burn <- funcyCtrlMbc@maxit
    CCO@eps <- funcyCtrlMbc@eps
    CCO@init <- init
    if(is.null(funcyCtrlMbc@seed))
        CCO@seed <- sample.int(10^9, size=1)
    else
        CCO@seed <- funcyCtrlMbc@seed
    CCO@loglikplot <- plotLoglik

    ptm <- proc.time()
    CCR <- getFCMM(CCDred,CCO)
    sysTime <- proc.time()-ptm

    ##funcyOut
    out <- new("funcyOutMbc")
    out@methodName <- "waveclust"
    out@kOut <- k
    out@time <- regTime
    out@dimBaseOut <- as.numeric(NA)
    cluster <- apply(CCR["Tau"], 1, which.max)
    names(cluster) <- NULL
    out@cluster<- cluster
    out@centers<-t(do.call(rbind,getwr.mu(CCR,CCO,CCDred)))
    out@props<-round(CCR@prop,4)
    out@dist2centers <- dist2centers(data, out@centers)
    out@cldist=makeClMat(out@dist2centers)
    out@calcTime <- sysTime 
    plotParams <- as.list(NA)
    ##funcyOutMbc
    #out@groupDimBase <- NA
    out@probs <- CCR@Tau
    out@prms <- list(Beta=CCR@Beta, Alpha=CCR@Alpha)
    aicbic <- getAICBIC(CCR, CCDred)
    out@AIC <- aicbic$AIC
    out@BIC <-aicbic$BIC
    out@logLik <- CCR@loglik
    out@nrIter <- as.numeric(NA)
       
    return(out)

}
