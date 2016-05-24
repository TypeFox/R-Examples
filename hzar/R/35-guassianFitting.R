

g.LLfuncA <- function(obsData, reqExp, muExp, varExp, tFunc, tArgs, tFixed,
           LLrejectedModel = -1e+08){
    baseFunc <- function(theta) 0;
    
    gLL <-
      guassianThetaLLExpF(distance=obsData$frame$dist,
                          sampleMean=obsData$frame$mu,
                          sampleVar=obsData$frame$var,
                          nEff=obsData$frame$nEff,
                          muExp=muExp,
                          varExp=varExp)
    gLL <-  step1VectorExpF(reqExp,gLL,LLrejectedModel)
    gLL <- eval(substitute(substitute(LLfunc,tF),list(LLfunc=gLL,tF=tFixed)))
    lapply(1:length(tArgs),function(x) bquote(theta[.(x)]))->tMap
    names(tMap) <- tArgs;
    body(baseFunc) <- simplify.exp(eval(substitute(substitute(LLfunc,eL),
                                      list(LLfunc=gLL,eL=tMap))))
    environment(baseFunc) <- .GlobalEnv
    baseFunc
}
g.LLfunc <- function(obsData, model,modelParam=splitParameters(model),
           tFixed=modelParam$fixed,
           tArgs=modelParam$init,
           LLrejectedModel = -1e+08){
    baseFunc <- function(theta) 0;
    model.req=model$req;
    model.gen=model$mFunc
    model.obsData=obsData
    param.fixed=tFixed
    frame=obsData$frame
    new.formals=c(formals(model.req)[names(tArgs)],tFixed)
    
    gLL <-
      guassianThetaLLExpF(distance=quote(frame$dist),
                          sampleMean=quote(frame$mu),
                          sampleVar=quote(frame$var),
                          nEff=quote(frame$nEff),
                          muExp=model$mu,
                          varExp=model$var)
    
    
    LLrejectedModel <- bquote(.(LLrejectedModel)+
                              .(req.edge.exp(body(model$req)[[2]])))
    gLL <-  step1VectorExpF(body(model$req)[[2]],gLL,LLrejectedModel)
    ## gLL <- eval(substitute(substitute(LLfunc,tF),list(LLfunc=gLL,tF=tFixed)))
    gLL <- ll.compile.theta(tArgs,tFixed,gLL)
    body(baseFunc) <-as.call(c(as.name("{"),
                               bquote(res <- .(gLL)),
                               expression(if(any(is.na(res))) print(theta)),
                               bquote(ifelse(is.na(res),.(LLrejectedModel),res))));
    llFunc=baseFunc
    formals(model.req) <- new.formals
    formals(model.gen) <- new.formals
    
    ## eval(substitute(substitute(LLfunc,eL),
    ##                                 list(LLfunc=gLL,eL=tMap)))
    ## body(baseFunc) <- substitute(evalq(LLfunc,envir=theta),list(LLfunc=gLL))
    ## environment(baseFunc) <- .GlobalEnv 
    baseFunc
}
ll.compile.theta <- function(tArgs,tFixed,gLL){

    lapply(1:length(tArgs),function(x) bquote(theta[[.(x)]]))->tMap
    names(tMap) <- names(tArgs);
    return(eval(substitute(substitute(LLfunc,tF),
                           list(LLfunc=gLL,tF=c(tFixed,tMap)))))
}


appThetaWalker <- function(theta,zF,tLower,tUpper,random=100,passCenter=FALSE){
## cat(as.numeric(theta))
tK <- names(theta)
tLower <- tLower[tK]
tUpper <- tUpper[tK]
##str(tLower)
##str(tUpper)
  tDelta <- lapply(tK,function(x)theta[[x]]-tLower[[x]])
names(tDelta) <- tK
##str(tDelta)
  n <- length(theta)
  z0 <- zF(theta)
  zM <- matrix(0,nrow=n,ncol=2)
  ##eZCa <- FALSE
  dimnames(zM) <- list(tK,c("A","C"))

eZC <- function(i,j){
    junk <- theta;
    if(j==1) junk[[i]]=theta[[i]]-tDelta[[i]]
    if(j==2) junk[[i]]=theta[[i]]+tDelta[[i]]
    zM[i,j]<<-zF(junk)
  }
  eZM <- function(){
    z0<<-zF(theta)
    for(iter in 1:n){
      eZC(iter,1);eZC(iter,2)
    }
  }
tCnt=as.list(rep(-1,n))
names(tCnt) <- tK
  rPmA <- function(i,d=0.5){
    ##cat(i,"r ",sep="")
    tDelta[[i]]<<-d*tDelta[[i]]
    tCnt[[i]] <<- tCnt[[i]]+1;
    eZC(i,1); eZC(i,2);
  }
  sPmU <- function(i,d=tDelta[[i]]) {
    ##cat(i,"> ",sep="")
    if(theta[[i]]>=tUpper[[i]]) stop("Moving outside appeture U")
    if(theta[[i]]+tDelta[[i]]>=tUpper[[i]])
      rPmA(i,0.25);
    if(theta[[i]]+2*tDelta[[i]]>=tUpper[[i]])
      rPmA(i,0.5);
    theta[[i]]<<-theta[[i]]+tDelta[[i]]
    zM[i,1]<<-z0
    z0<<-zM[i,2]
    eZC(i,2)
  }
  sPmD <- function(i,d=tDelta[[i]]) {
    ##cat(i,"< ",sep="")
    if(theta[[i]]<=tLower[[i]]) stop("Moving outside appeture D")
    if(theta[[i]]-tDelta[[i]]<=tLower[[i]])
      rPmA(i,0.25);
    if(theta[[i]]-2*tDelta[[i]]<=tLower[[i]])
      rPmA(i,0.5);
    theta[[i]]<<-theta[[i]]-tDelta[[i]]
    zM[i,2]<<-z0
    z0<<-zM[i,1]
    eZC(i,1)
  }
  for(iter in tK){
    junk <- tUpper[[iter]]-theta[[iter]]
    if(tDelta[[iter]]==0){
      tDelta[[iter]]= junk
      rPmA(iter,0.1)
      sPmU(iter)
    } else if(junk != 0) {
      tDelta[[iter]]=min(junk,tDelta[[iter]])
    } else {
      rPmA(iter,0.1)
      sPmD(iter)
    }
    rPmA(iter,0.1)
  }
  eZM();
  wUpL=list()
wAL=as.list(rep(FALSE,n))
names(wAL) <- tK
  wU <- function(i){
    if(!wAL[[i]]) {
      wUpL[[i]]<<-TRUE
      wAL[[i]]<<-TRUE
    }else if(!wUpL[[i]]) wR(i)
  }
  wD <- function(i){
    if(!wAL[[i]]) {
      wUpL[[i]]<<-FALSE
      wAL[[i]]<<-TRUE
    }else if(wUpL[[i]]) wR(i)
  }
  wR <- function(i){
    wUpL[[i]]<<-!wUpL[[i]]
    rPmA(i) 
  }
  ## zHill <- function(a,b,c) (2*zF(b)) > (zF(a)+zF(c))
  ## zTop <- function(a,b,c) (zF(b)>zF(a))&(zF(b)>zF(c))
  zMCk <- function()
    all(2*z0 > zM[,1] +zM[,2])
  zMCtop <- function()
    (z0 > zM[,1])&(z0>zM[,2])
  zMMd <- function()
    tK[(2*z0 <= zM[,1] +zM[,2]) &(zM[,1]!=zM[,2])]
  
  zPmW <- function(i){
    ## zA <- zM[i,1];
##     zB <- z0
##     zC <- zM[i,2];
    while(!(2*z0> zM[i,1]+zM[i,2])){
      #Move up
      if( zM[i,1]> zM[i,2]){
        sPmD(i)
        wD(i)
      }else if( zM[i,2]> zM[i,1]){
        sPmU(i)
        wU(i)
      }else{
        #Flat parameter space. Off edge?
        eZM();
        if((theta[[i]]-tDelta[[i]]*2>tLower[[i]])&&
           (theta[[i]]+tDelta[[i]]*2<tUpper[[i]]))
          rPmA(i,2)
        else
          rPmA(i,min(theta[[i]]-tLower[[i]],
                     tUpper[[i]]-theta[[i]])*0.9/tDelta[[i]])
        ## lapply(tK[tK!=i],rPmA)
        return(-1)
      }
    }
    if(tCnt[[i]]>0) return(0)
    ## Curve good
    if(zM[i,2]>z0){
      sPmD(i)                         #Need to go down to inflection
      wD(i)
    } else if(zM[i,1] > z0){
      sPmU(i)                         #Need to go up to inflection
      wU(i)
    } else {
      rPmA(i,0.9)                           #At crown, shrink appeture
      eZM();
      return(1)
    }
    eZM();
    return(0)
  }
## print(z0)
##   print(cbind(t=theta,
##               d=tDelta))
## print(zM)
  for(zzz in tK){
    sapply(zMMd(),zPmW)
    if(zMCk()) break;
  }
## print(z0)
##   print(cbind(t=theta,
##               d=tDelta))
## print(zM)

  covTest <- function(){
    hM <- NULL;
    cV <- NULL;
    test <- NULL;
    try({
      hM <- -appHessian(theta,zF,tDelta)
      if(all(diag(hM)>0)||all(diag(hM)<0)){
        cV <- solve(hM);
      
    ##    cat(diag(cV)>0,"\n")
        if(all(diag(cV)>0))
          test <- chol(cV);
      }
    }, silent=TRUE)
    !is.null(test)
  }
        
  for(zzz in tK){
    
    if(all(zMCtop()))break;
    if(covTest())break;
    sapply(tK,zPmW)
    ##cat(" || ")
  }
##cat("\n")

## print(z0)
##   print(cbind(t=theta,
##               d=tDelta))
## print(zM)
## print(solve(-appHessian(theta,zF,tDelta)))
  if(passCenter)
  return(list(cov=solve(-appHessian(theta,zF,tDelta)),
              center=theta));
  
  return(solve(-appHessian(theta,zF,tDelta)));
  
}
appHScale <- function(cV,tLower,tUpper){
  for(iter in colnames(cV)){
    maxV <- (tUpper[[iter]]-tLower[[iter]])^2
    if(cV[iter,iter] > maxV){
      r <- sqrt(maxV/cV[iter,iter])/2
      for(jI in colnames(cV)){
        cV[jI,iter]=cV[jI,iter]*r
        cV[iter,jI]=cV[iter,jI]*r
      }
    }   
  }
  return(cV)
}
appThetaWalkerR <- function(theta,zF,tLower,tUpper,random=1000,passCenter=FALSE){
  while((random <- (random -1) )>0){
    try({
    res <- appThetaWalker(theta,zF,tLower,tUpper,10,passCenter)
    if(passCenter){
      if(all(diag(res$cov)>0)) {
        chol(res$cov);
        res$cov <- appHScale(res$cov,tLower,tUpper)
        return(res);}
    }else{
      if(all(diag(res)>0)) {
        chol(res);
        res <- appHScale(res,tLower,tUpper)
        return(res);}
    } }, silent=TRUE)
    for(iter in names(theta)){
      theta[[iter]] <- runif(1,tLower[[iter]],tUpper[[iter]])
    }
  }
}
appHessian <- function(theta,zF,tDelta=as.list(0.05*as.numeric(theta))){
  zFunc <- function(a,b,c,d)
    zF(a)+zF(d)-zF(b)-zF(c);
  n=length(theta)
  res <- matrix(0,nrow=n,ncol=n);
  for(i in 1:n){
    for(j in i:n){
      
      if(i==j){
        t0=t2=theta
        t0[[i]]=theta[[i]]-tDelta[[i]]
        t2[[i]]=theta[[i]]+tDelta[[i]]
        res[i,i]=zFunc(t0,theta,theta,t2)/(tDelta[[i]]^2);
    
      }else{
        t00=t02=t20=t22=theta
        t00[[i]]=theta[[i]]-tDelta[[i]]
        t02[[i]]=theta[[i]]-tDelta[[i]]
        t20[[i]]=theta[[i]]+tDelta[[i]]
        t22[[i]]=theta[[i]]+tDelta[[i]]
        t00[[i]]=theta[[j]]-tDelta[[j]]
        t20[[i]]=theta[[j]]-tDelta[[j]]
        t02[[i]]=theta[[j]]+tDelta[[j]]
        t22[[i]]=theta[[j]]+tDelta[[j]]
        res[j,i]=res[i,j]=
          zFunc(t00,t20,t02,t22)/(tDelta[[i]]*tDelta[[j]]);
      }
    }
  }
  rownames(res) <- names(theta)
  colnames(res) <- names(theta)
  
  
  return(res);

}
naiveHessian <- function(theta,zF,k=0.05){
  zFunc <- function(a,b,c,d)
    zF(a)+zF(d)-zF(b)-zF(c);
  n=length(theta)
  res <- matrix(0,nrow=n,ncol=n);
  for(i in 1:n){
    for(j in i:n){
      
      if(i==j){
        t0=t2=theta
        t0[[i]]=theta[[i]]*(1-k)
        t2[[i]]=theta[[i]]*(1+k)
        res[i,i]=zFunc(t0,theta,theta,t2)/(k^2*theta[[i]]^2);
    
      }else{
        t00=t02=t20=t22=theta
        t00[[i]]=theta[[i]]*(1-k)
        t02[[i]]=theta[[i]]*(1-k)
        t20[[i]]=theta[[i]]*(1+k)
        t22[[i]]=theta[[i]]*(1+k)
        t00[[i]]=theta[[j]]*(1-k)
        t20[[i]]=theta[[j]]*(1-k)
        t02[[i]]=theta[[j]]*(1+k)
        t22[[i]]=theta[[j]]*(1+k)
        res[j,i]=res[i,j]=
          zFunc(t00,t20,t02,t22)/(k^2*theta[[i]]*theta[[j]]);
      }
    }
  }
  rownames(res) <- names(theta)
  colnames(res) <- names(theta)
  
  
  return(res);

}

hzar.first.fitRequest.gC <- function(gModel,obsData,verbose=TRUE){
  if (verbose) {
    mcmcParam <- cfg.hzar.default.mcmc
  } else {
    mcmcParam <- cfg.hzar.quiet.mcmc
  }
  modelParam <- splitParameters(gModel$parameterTypes);
  
  LLfunc <- g.LLfunc(obsData,gModel,modelParam)

##   Figure out a better covariance matrix method.

##   I currently first try an approximate Hessian around the initial
##   parameters.  If that fails, I then try to directly approximate the
##   covariance matrix using hzar.cov.rect.  If that fails to work, I
##   should walk back over the parameter space (possibly using the values
##   pulled by hzar.cov.rect) starting at the maximum likelihood,
##   searching for the first set of parameters that the first method
##   produces useful values.

  cV <- NULL;
  junk <- NULL
  altInit <- modelParam$init
  if(verbose) cat("a")
  try({ junk <-  solve(-naiveHessian(modelParam$init,LLfunc));
        junk <- appHScale(junk,modelParam$lower, modelParam$upper)
      },silent=TRUE)
  if(!is.null(junk)){
    
      try(if(all(diag(junk)>0)){chol(junk); cV <- junk;},silent=TRUE)
   
  }
  else
    if(verbose) cat("A")
  if(is.null(cV)){
    if(verbose) cat("b")
    try({
      junk <- appThetaWalkerR(modelParam$init,
                             LLfunc, modelParam$lower, 
                            modelParam$upper, random = 1000,
                            passCenter=TRUE);
      cV <- junk$cov
      modelParam$init  <- junk$center
      
    })
    junk <- NULL
    ## print(cV)
    try(junk <- chol(cV),silent=TRUE)
    if(is.null(junk)){
      if(verbose) cat("c")
      modelParam$init <- altInit
      try(  junk <-  solve(-naiveHessian(modelParam$init,LLfunc)),silent=TRUE)
      if(!is.null(junk))
        try({ chol(junk); cV <- junk; },silent=TRUE)
      else
        if(verbose) cat("C")
    ## data.mat<-list(dTheta=prod(abs(as.numeric(param.upper)
##                      -as.numeric(param.lower)))
      ##                    / random,
##                    data=fitter.gen.rParam.uniform(param.lower,
##                      param.upper,
##                      random));
    }
    
  }
  
  cat("\n")
  return(hzar.make.fitRequest(modelParam, cV, LLfunc, 
                              mcmcParam))
}
