
thEvalFExp <- function( target, exp, val)
  bquote(.(target) <- .(eval(substitute(substitute(a,b),
                                        list(a=exp,b=val)))))
optMuPExpF <- function( sampleMean, nEff,distance,
                       fExp, varExp,
                       tMuL=quote(muL),
                       tMuR=quote(muR)){
  fTa <- quote(fXA)
  fTd <- quote(fXD)
  varT <- quote(vX)
  cLT <- quote(cL)
  cRT <- quote(cR)
  aLT <- quote(aL)
  aRT <- quote(aR)
  aMT <- quote(aM)
  tDen <- quote(aLaRaMM)
  
  
  c(thEvalFExp(fTa,fExp,list(x=distance)),
    bquote(.(fTd) <- 1 - .(fTa)),
    thEvalFExp(varT,varExp,list(x=distance)),
    bquote(.(cLT) <- sum( .(sampleMean)*.(nEff)*.(fTd)/ .(varT))),
    bquote(.(cRT) <- sum( .(sampleMean)*.(nEff)*.(fTa)/ .(varT))),
    bquote(.(aLT) <- sum( .(nEff)*.(fTd)^2/ .(varT))),
    bquote(.(aRT) <- sum( .(nEff)*.(fTa)^2/ .(varT))),
    bquote(.(aMT) <- sum( .(nEff)*.(fTd)*.(fTa)/ .(varT))),
    bquote(.(tDen) <- .(aRT)*.(aLT)-.(aMT)^2),
    bquote(.(tMuL) <- (.(cLT)*.(aRT)-.(cRT)*.(aMT))/.(tDen)),
    bquote(.(tMuR) <- (.(cRT)*.(aLT)-.(cLT)*.(aMT))/.(tDen))
    )
}

mu2fExp <- function(muExp, tMuL=quote(muL), tMuR=quote(muR))
  simplify.exp(eval(substitute(substitute(a,b),list(a=muExp,b=list(muL=0,muR=1)))))

optMuPExpF.fit <- function(fitRequest,
                           model=cline.extract.modelPrep(fitRequest)$model,...){
  if(!all(c("mu","var") %in% names(model)))
    return(list());
  fixed=cline.extract.modelPrep(fitRequest)$modelParam$fixed
  
  optMuPExpF(distance=quote(frame$dist),
             sampleMean=quote(frame$mu),
             nEff=quote(frame$nEff),
             fExp=mu2fExp(eval(substitute(substitute(a,fixed),
               list(a=model$mu)))),
             varExp=simplify.exp(eval(substitute(substitute(a,fixed),
               list(a=model$var)))),
             ...)
}

optimizeTraceF <- function(fitRequest){
  resF <- function(mcmcRaw) {as.data.frame(mcmcRaw)}
  params <- names(cline.extract.modelPrep(fitRequest)$modelParam$init)
  params <- params[!(params %in% c("muL","muR"))]
  repTheta <- lapply(params,function(x) bquote(data.param[[.(x)]][iter]))
  names(repTheta) <- params
  ## str(repTheta)
  optMuPExpF.fit(fitRequest,
                 tMuL=quote(data.param$muL[iter]),
                 tMuR=quote(data.param$muR[iter]))->optExp
  optExp <- lapply(optExp,
                   function(x) eval(substitute(substitute(a,b),
                            list(a=x,b=repTheta))))
  as.call(c(as.name("{"),
            c(quote(data.param<-as.data.frame(mcmcRaw)),
              bquote(for(iter in 1:nrow(data.param))
                     .(as.call(c(as.name("{"),optExp)))),
              quote(data.param)) ) ) -> body(resF)
  environment(resF) <- list2env(list(frame=hzar.extract.obsData(fitRequest)$frame),
                                parent=.GlobalEnv)
  resF
}

optimizeTrace <- function(fitRequest,...){
  ##print(
    traceF <- optimizeTraceF(fitRequest)
  ##  )
  
  if(inherits(fitRequest,"hzar.fitRequest")){
    if(!all(c("mcmcRaw") %in% names(fitRequest)))
      stop("No trace to optimize")
    
    tempPar <- attr(fitRequest$mcmcRaw, "mcpar")
    fitRequest$mcmcRaw <- as.mcmc(traceF(fitRequest$mcmcRaw))
    tempPar -> attr(fitRequest$mcmcRaw, "mcpar")
    return(fitRequest)
  } else if (inherits(fitRequest,"hzar.dataGroup")){
    tempPar <- attr(fitRequest$data.mcmc, "mcpar")
    
    fitRequest$data.param <-traceF(fitRequest$data.mcmc)
    fitRequest$data.mcmc <- as.mcmc(fitRequest$data.param )
    tempPar -> attr(fitRequest$data.mcmc, "mcpar")
    return(hzar.make.dataGroup(data.mcmc=fitRequest$data.mcmc,
                               llFunc=fitRequest$llFunc,
                               data.param=fitRequest$data.param,
                               obsData=fitRequest$obsData,
                               ...))
  }
}



profile.model <- function(parent.model, parameter,fixed.value){
  if(length(fixed.value)>1)
    return(lapply(fixed.value,profile.model,parent.model=parent.model,parameter=parameter))
  res=parent.model
  meta.init(res)[[parameter]] <- fixed.value
  meta.fix(res)[[parameter]] <- TRUE
  meta.tune(res) <- 1.1
  res
}


profile.fitRequest <- function(model,obsData, old.mcmc, mcmcParam =hzar.make.mcmcParam(1e5,1e4,5e4,100) ){
 ##print(summary(old.mcmc)) 
  mdlParam<-splitParameters(model$parameterTypes);
  if(inherits(obsData,"guassSampleData1D")){
    clineLLfunc <- g.LLfunc(obsData, model,modelParam=mdlParam)
  }else if(inherits(obsData,"clineSampleData1D")){
    clineLLfunc <- freq.LLfunc(obsData,model,mdlParam$init,mdlParam$fixed)
  } else {
    stop(paste("Failed to process obsData class: <",paste(class(obsData),collapse=", "),">"))
  }
  covMatrix<-NULL;
  if(all(c("model.LL")%in%colnames(old.mcmc) ))
    mcmc.LL <- old.mcmc[,"model.LL"]
  else
    mcmc.LL <- NULL
  if((is.numeric(old.mcmc)||is.data.frame(old.mcmc))&&nrow(old.mcmc)>20){
    old.mcmc <- old.mcmc[,names(mdlParam$init),drop=FALSE]
    nGen <- dim(old.mcmc)[[1]]
    if(is.null(mcmc.LL)){
      nSam <- min(nGen,5e3)
      sCut <- sample(nGen,nSam)
      mcmcSubset<-old.mcmc[sCut, ,drop=FALSE];
      subLL<-hzar.eval.clineLL(mcmcSubset,clineLLfunc);
    }else{
      mcmcSubset <- old.mcmc
      subLL <- mcmc.LL
    }
        covData <- NULL
    if(sum(unique(subLL)>max(subLL-4))<0.95*nGen){
      try({
        if(sum(unique(subLL)>max(subLL-4))>1e3){
          wt <- subLL[subLL>max(subLL-4),]
          wt <- exp(wt-max(wt))
          covData <- cov.wt(mcmcSubset[subLL>max(subLL-4), ,drop=FALSE],wt)
        }else if(sum(unique(subLL)>max(subLL-4))>1e2){
          covData<-hzar.cov.mcmc(clineLLfunc,mcmcSubset[subLL>max(subLL-4), ,drop=FALSE],passCenter=TRUE);
        }else if(sum(unique(subLL)>max(subLL-10))>1e2){
          covData<-hzar.cov.mcmc(clineLLfunc,mcmcSubset[subLL>max(subLL-10), ,drop=FALSE],passCenter=TRUE);
        }else {
          covData<-hzar.cov.mcmc(clineLLfunc,mcmcSubset,passCenter=TRUE);
        }
      },silent=TRUE)
    }
    if(is.null(covData))
      try(covData <- cov.wt(mcmcSubset[subLL>max(subLL-4), ,drop=FALSE]))
    if(!is.null(covData)){
      covMatrix<-covData$cov;
      new.center<-covData$center[names(mdlParam$init)];
      if(clineLLfunc(new.center)>-1e6)
        mdlParam$init <- new.center;
    }
    if(!is.null(covMatrix)){
      junk <- covMatrix;
      covMatrix <- NULL;
      try(if(all(diag(junk)>0)){
        chol(junk);
        covMatrix <- junk;
      },silent=TRUE)
    }
    if(is.null(covMatrix)){
      try({ junk <-  solve(-naiveHessian(mdlParam$init,clineLLfunc));
            junk <- appHScale(junk,mdlParam$lower, mdlParam$upper);
            if(all(diag(junk)>0)){
              chol(junk);
              covMatrix <- junk;}
          },silent=TRUE)
    }
    if(is.null(covMatrix)){
      try({
        junk <- appThetaWalkerR(mdlParam$init,
                                clineLLfunc,
                                mdlParam$lower, 
                                mdlParam$upper, random = 1000,
                                passCenter=TRUE);
        covMatrix <- junk$cov
        mdlParam$init  <- junk$center
        
      })
    }
  }else{
    try(  covMatrix<-hzar.cov.rect(clineLLfunc,mdlParam$lower,mdlParam$upper,random=1e4));
  }
  return(hzar.make.fitRequest(mdlParam,covMatrix,clineLLfunc,mcmcParam));
}

hzar.profile.dataGroup <- function(dG, parameter,
                                   pVals=NULL,pDivs=NULL, nDiv=20,appeture=NULL,
                                   doPar=FALSE,...){
  dG.trace <- as.data.frame(hzar.mcmc.bindLL(dG))

  ##Get parameter space
  param.trace <- dG.trace[[parameter]]
  ##Get parameter intervals
  if(is.list(pDivs)){
    param.divs <- pDivs
    nDiv <- length(pDivs)
  }else if(is.numeric(pDivs)){
    param.divs <- lapply(2:length(pDivs),function(x) pDivs[c(-1,0)+x])
    nDiv <- length(param.divs)
  }else if(is.numeric(pVals)){
    pVals=sort(unique(pVals))
    if(!all(is.numeric(appeture),length(appeture)==1,appeture>0))
      appeture=median(pVals[2:length(pVals)]-pVals[1:(length(pVals)-1)])
    param.divs <- lapply(pVals,function(x) c(x-appeture,x+appeture))
    nDiv <- length(param.divs)
  }else{
    param.divs <- levels(equal.count(param.trace,nDiv))
  }
  param.subsets <- lapply(param.divs,function(x) which(param.trace>=x[1] &param.trace<=x[2]))
  if(is.numeric(pVals)){
    if(length(pVals)!=nDiv) stop("length mismatch between pVals and pDivs")
  } else{
    pVals=lapply(param.divs,mean)
  }
  
  
  ##Get parent model
  parent.model <- cline.extract.modelPrep(dG)$model
  
  obsData <- hzar.extract.obsData(dG)

  ##Build multi fit request

  profileFitR <- list();

  lP <- foreach(param.divs=param.divs,param.val=pVals,
                mcmc.subset=lapply(param.subsets,
                  function(x) dG.trace[x, ,drop=FALSE]),
                .packages="hzar")
  exec <- expression({
    temp.model <- profile.model(parent.model,parameter,param.val)
    ##Get non-profile parameters
    param.fixed <- as.logical(meta.fix(temp.model))
    param.names <- names(param.fixed)[!param.fixed]
 
    param.init <- mcmc.subset[which.max(mcmc.subset$model.LL),param.names,drop=FALSE]
    for(tk in param.names)meta.init(temp.model)[[tk]] <- param.init[1,tk]

    profile.fitRequest(temp.model,obsData,mcmc.subset)
  })
  ## deps=c("
  if(doPar)
    profileFitR <- lP %dopar% eval(exec)
  else
    profileFitR <- lP %do% eval(exec)
  hzar.multiFitRequest(profileFitR,...)
}
