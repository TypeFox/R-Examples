## This is the new way of fitting the cline model, in an attempt to
## streamline the process. I want to break the code into smaller
## chunks by taking advantage of internal methods.  That way I can
## provide complex behavior to the end user with minimal code.
## Although technically the code will be more complex in structure,
## each method should have a clear meaning.

## In writing this, I plan on starting with biggest atomic slice of
## the fitting process, the call to MCMCmetrop1R.

hzar.doFit <- function(fitRequest){
  result<-NULL;                        
  useParam=fitRequest$mcmcParam;
  mdlParam=fitRequest$modelParam;
  try({
    result<-
      MCMCmetrop1R(fun=fitRequest$llFunc, logfun="true", force.samp=TRUE,
                   mcmc=useParam$chainLength, burnin=useParam$burnin,
                   verbose=useParam$verbosity, thin=useParam$thin,
                   theta.init=mdlParam$init, tune=as.numeric(mdlParam$tune),
                   V=fitRequest$cM, seed=useParam$seed,
                   optim.control=list(fnscale=-1,trace=0,REPORT=10,maxit=5000),
                   optim.method="L-BFGS-B",
                   optim.lower=useParam$lower,
                   optim.upper=useParam$upper );
    
    colnames(result)<-names(mdlParam$init);
  });
  fitRequest$mcmcRaw<-result;
  attr(fitRequest,"fit.run")<-TRUE;
  attr(fitRequest,"fit.success")<-FALSE;
  if(identical(is.null(result),TRUE)){
    warning("Fitting failed.");
  } else {
    attr(fitRequest,"fit.success")<-TRUE;
  }
  return(fitRequest);
}

## fitRequest has the following base structure:
hzar.make.fitRequest <-
  function(modelParameters,             #output of splitParameters
           covMatrix,                   #covariance matrix to use
           clineLLfunc,                 #method to get LL given theta
           mcmcParameters,              #mcmc attributes & controls
           mcmcRaw=NULL,                #what MCMCmetrop1R returns
           fit.run=FALSE,               #has this request been run?
           fit.success=FALSE            #has the run succeeded?
           ){
    fitRequest<-list(modelParam=modelParameters,cM=covMatrix,
                     llFunc=clineLLfunc,mcmcParam=mcmcParameters,
                     mcmcRaw=mcmcRaw);
    class(fitRequest)<-"hzar.fitRequest";
    attr(fitRequest,"fit.run")<-fit.run;
    attr(fitRequest,"fit.success")<-fit.success;
    return(fitRequest);
  }

## I already have splitParameters to make modelParameters.

## mcmcParameters should have a simple structure:
##
## chainLength, burnin, verbosity, thin, seed(seedStreamChannel,
## useSeedStream, mersenneSeed, lecuyerSeed)
hzar.make.mcmcParam <-
  function(chainLength, burnin, verbosity, thin,
           seedStreamChannel=1, useSeedStream=TRUE,
           mersenneSeed=12345,
           lecuyerSeed=rep(12345,6)){
    mcmcParam<-list(chainLength=chainLength,burnin=burnin,
                    verbosity=verbosity,thin=thin);
    mcmcParam$seed<-mersenneSeed;
    if(useSeedStream){
      mcmcParam$seed<-list(lecuyerSeed,seedStreamChannel)
    }
    return(mcmcParam);
  }



## I need a method to make clineLLfunc.  I will start with the
## original cline log likelihood function, wrapped with a basic method
## call.

## I need two llFuncs, one for Maximum Likelihood and one for Bayesian
## analysis. Also, I need to recognize that under the old process, the
## fixed parameters are magically inserted during engine startup to
## simplify the method calls.
hzar.make.clineLLfunc.old.ML <-
  function(param.free.names, param.fixed, #lists made by splitParam
           param.check.func,              #func to check parameters
           meta.cline.func,               #func to generate cline func
           model.LL,                      #returns LL of cline func
           LLrejectedModel=-1E8           #when rejecting, return this
           ){
    ## First, assign references locally.
    model.req<-param.check.func;
    model.gen<-meta.cline.func;
    eval.clineLL<- model.LL;
    obsData=as.list(environment(model.LL))$obj;
    myRejectionLL<-LLrejectedModel;
    ## Second, modify function signitures to account for fixed
    ## parameters.
    old.formals<-formals(model.gen);
    if(length(old.formals)!=(length(param.free.names)+
                             length(param.fixed))){
      warning("The length of the method formals does not match the length of the parameters supplied.");
    }
    ttt.formals<-old.formals[param.free.names];
    names(ttt.formals)<-param.free.names;
    new.formals<-c(ttt.formals,param.fixed);
    formals(model.req)<-new.formals;
    formals(model.gen)<-new.formals;
    
    ## references: theta, meta.model, obsData, myRejectionLL
    ##
    ## in-depth:  theta, meta.model$req, myRejectionLL, meta.model$func, 
    ## obsData$model.LL
    
    llFunc<-function(theta){
      if(! do.call(model.req,as.list(theta))) return(myRejectionLL);
      model=do.call(model.gen,as.list(theta));
      result<-eval.clineLL(model);
      if(identical(is.finite(result),TRUE))
        return(result);
      return(myRejectionLL);
    }
    
    return(llFunc);
  }

hzar.make.clineLLfunc.old.bayes <-
  function(param.free.names, param.fixed, #lists made by splitParam
           param.check.func,              #func to check parameters
           meta.cline.func,               #func to generate cline func
           model.LL,                      #returns LL of cline func
           prior.LL,                      #returns LL of priors
           LLrejectedModel=-1E8           #when rejecting, return this
           ){

    ## First, assign references locally.
    model.req<-param.check.func;
    model.gen<-meta.cline.func;
    eval.priorLL<-prior.LL;
    eval.clineLL<- model.LL;
    myRejectionLL<-LLrejectedModel;
    ## Second, modify function signitures to account for fixed
    ## parameters.
    old.formals<-formals(model.gen);
    if(length(old.formals)!=(length(param.free.names)+
                             length(param.fixed))){
      warning("The length of the method formals does not match the length of the parameters supplied.");
    }
    
    ttt.formals<-old.formals[param.free.names];
    names(ttt.formals)<-param.free.names;
    new.formals<-c(ttt.formals,param.fixed);
    formals(model.req)<-new.formals;
    formals(model.gen)<-new.formals;
    formals(eval.priorLL)<-new.formals;
    
    ## references: theta, meta.model, obsData, myRejectionLL
    ##
    ## in-depth: theta, meta.model$req, myRejectionLL, meta.model$func,
    ## meta.model$prior, obsData$model.LL
    llFunc<-function(theta){
      if(! do.call(model.req,as.list(theta))) return(myRejectionLL);
      model=do.call(model.gen,as.list(theta));
      thetaLL=do.call(eval.priorLL,as.list(theta));
      
      result<-eval.clineLL(model)+thetaLL;
      if(identical(is.finite(result),TRUE))
        return(result);
      return(myRejectionLL);
    }
    
    return(llFunc);
  }


## I need method(s) to generate covMatrix.

## I will need helper functions, specifically one to generate a vector
## of likelihood values given a data frame of parameter values. I want
## this to be capable of executing in parallel.

## require(foreach);

fitter.wedgeSlice <- function (count,slice.size=1000){
  res=list();
  
  if(count>slice.size)
    res<-lapply(1:((count-1)%/%slice.size)-1,
                function(x,y=slice.size) { 1:y + x*y });
  return(c(res,list(1:(1+(count-1)%%slice.size)
                    +slice.size*((count-1)%/%slice.size))));
}

hzar.eval.clineLL <- function(data, llFunc,doPar=FALSE){
  ## print("A");
  slices<-fitter.wedgeSlice(dim(data)[[1]]);
  ## cat("Eval wedge size:",object.size(slices),"\n");
  tttIndex=NULL;
  useFunc=llFunc;
  useData=data;
  if(doPar){
    
    result<-foreach(tttIndex=slices,
                    .combine=rbind) %dopar% {
                      data.frame(model.LL=as.numeric(lapply(tttIndex,function(x)
                                   {useFunc(useData[x,]); })));}
  }else{
    
    result<-foreach(tttIndex=slices,
                    .combine=rbind) %do% {
                      data.frame(model.LL=as.numeric(lapply(tttIndex,function(x)
                                   {useFunc(useData[x,]); })));}
  }
  ## cat("Eval result size:",object.size(result),"\n");
  return(result);
}
           
fitter.gen.rParam.uniform<-function(param.lower,param.upper,count=1000){
  low=NULL; high=NULL;
  raw<-foreach(low=param.lower,
               high=param.upper,
               .combine=cbind) %do% {runif(count,low,high)};
  result<-as.data.frame(raw);
  colnames(result)<-names(param.lower);
  names(result)<-names(param.lower);
  return(result);
}

## So, what do I need to know to generate the matrix?
## Work backwards?

## cov.wt requires a data frame of sampled parameters, weights
## (Likelihood, intregration differentials, scaling), and possibly a
## pre-calculated center.

## data    := sampled parameters

## data.wt := weights needs: (data, clineLLfunc, data.dP

## data.LL := Likelihood; needs (clineLLfunc)
## data.LL =  hzar.eval.clineLL(data, clineLLfunc)
## data.dP := area given to each sampled point
## scaling := 1/(sum(exp(data.LL)*data.dP))

fitter.getCovWeights <- function(data,clineLLfunc, data.dP){
  data.LL <- hzar.eval.clineLL(data,clineLLfunc);
  return(exp(data.LL)*data.dP/sum(exp(data.LL)*data.dP));
}

## I can assume I have clineLLfunc, so I just need data and data.dP. I
## could just sample a rectangular area.

fitter.gen.samples.rect <- function(param.lower, param.upper, pDiv=11){
  param.names<-names(param.lower);
  nParam<-length(param.names);
  deltas=(as.numeric(param.upper[param.names])-
          as.numeric(param.lower[param.names]))/(pDiv-1);
 ##  print(deltas);print(param.names);
  grid.formals<-matrix(nrow=pDiv,ncol=nParam,
                       data=rep(deltas,each=pDiv))* rep(0:(pDiv-1),nParam)+
                         rep(as.numeric(param.lower[param.names]),each=pDiv);
  colnames(grid.formals)<-param.names;
  ## names(grid.formals)<-param.names;
  ## print(as.list(as.data.frame(grid.formals)));
  result<-list(dTheta=prod(abs(deltas)));
  result$data<-do.call(expand.grid,as.list(as.data.frame(grid.formals)));
  ## print(class(result$data))
  ## grid.formals<-names(param.lower
  return(result);
}

hzar.cov.rect<-function(clineLLfunc,param.lower,param.upper,pDiv=11,random=0,passCenter=FALSE){
  ## print("A");
  cat("0")
  ## Check to make sure we aren't about to do something stupid.
  ## Yes, that means that models with over ten million parameters
  ## will fail spectacularly. We hope that you would have
  ## considered writing your own software for problems of such
  ## scale.
  if(random>1e5){
    stop("Covariance matrix calculation with random sampling requested with far too many samples.  Stopping.")
  }
  if(random>0){
    cat("a")
    data.mat<-list(dTheta=prod(abs(as.numeric(param.upper)
                     -as.numeric(param.lower)))
                   / random,
                   data=fitter.gen.rParam.uniform(param.lower,
                     param.upper,
                     random));
  }else{
    cat("b")
    ## Stupidity check.  See above
    if((pDiv^length(param.lower))>1e6){
      ## This is recoverable by switching to random sampling.
      warning("Covariance matrix calculation requested for too complex of a lattice structure. Switching to random sampling using ten thousand samples.");
      return(hzar.cov.rect(clineLLfunc,
                           param.lower,
                           param.upper,
                           random=1e4,
                           passCenter=passCenter));
    }
    cat("B")
    data.mat<-fitter.gen.samples.rect(param.lower,param.upper,pDiv);
  }
  param.names<-names(data.mat$data);##print(names(data.mat));
  ## print("A");
  ##data.wt<-fitter.getCovWeights(data.mat$data,clineLLfunc,data.mat$dTheta);
  cat("c",sprintf("%0.1e",data.mat$dTheta))
  data.wt<-hzar.eval.clineLL(data.mat$data,clineLLfunc);
  data.mat$data<-data.mat$data[data.wt>-1e7, ,drop=FALSE];
  data.wt<-data.wt[data.wt>-1e7];
  cat("d")
  MIN.DATA<-(1+length(param.upper))*4
  if(length(data.wt)<MIN.DATA){
    ## need more samples... recurse.
    
    if(length(data.wt)>9){
      cat("-")
      param.A <- apply(data.mat$data,2,extendrange)
      
      if(!any(param.A[1,]==param.A[2,])){
        param.lower <- param.A[1,]
        param.upper <- param.A[2,]
        random=random/2
      }
    }
    
    if(random>0){
      cat("Dr")
      ## Double the amount of sampling
      return(hzar.cov.rect(clineLLfunc,
                           param.lower,
                           param.upper,
                           random=2*random,
                           passCenter=passCenter));
    }else{
      cat("Dp")
      ##Increase the resolution slightly.
      return(hzar.cov.rect(clineLLfunc,
                           param.lower,
                           param.upper,
                           pDiv=pDiv+1,
                           passCenter=passCenter));
    }
  }
  while(sum(data.wt>-723)<MIN.DATA){
    cat("+")
    ## insufficient data in finite range
    if(sum(data.wt>609)>0){
      ## scaling won't fix the problem, so more samples needed.
      ## Recurse.
      
      if(sum(data.wt>-1000)>9){
        cat("-")
        data.mat$data <- data.mat$data[data.wt>-1000, ,drop=FALSE]
        data.wt <- data.wt[data.wt>-1000]
      }
      
      cat("-")
      param.A <- apply(data.mat$data,2,extendrange)
      
      if(!any(param.A[1,]==param.A[2,])){
        param.lower <- param.A[1,]
        param.upper <- param.A[2,]
      }
      if(random>0){ 
        cat("Er")
        ## Double the amount of sampling
        return(hzar.cov.rect(clineLLfunc,
                             param.lower,
                             param.upper,
                             random=random*2,
                             passCenter=passCenter));
      }else{
        
        cat("Ep")
        ##Increase the resolution slightly.
        return(hzar.cov.rect(clineLLfunc,
                             param.lower,
                             param.upper,
                             pDiv=pDiv,
                             passCenter=passCenter));
      }
    }
    ## iteratively shift likelihood space to bring samples into finite range.
    if(sum(data.wt>-723)==0){
      data.wt<-data.wt-max(data.wt);
    }else {data.wt<-data.wt+100}
  }
  ## print("A");

  cat("f")
  VDATA<-cov.wt(x=cbind(data.mat$data,model.LL=data.wt),wt=exp(data.wt))
  ##   *data.mat$dTheta)
  VMATRIX<-VDATA$cov;
  ## diag(1/(VMATRIX["model.LL",]))->counter.inv;
  ## diag(sign(VMATRIX["model.LL",]))->counter.inv;
  ## dimnames(counter.inv)<-list(rownames(VMATRIX),colnames(VMATRIX));
  ## counter.inv2<-counter.inv/sqrt(counter.inv["model.LL","model.LL"]);
  ## mat.scaled<-counter.inv%*%VMATRIX%*%counter.inv;
  mat.scaled<-VMATRIX;
  
  cat("g")
  if(passCenter)
    return(list(cov=mat.scaled[param.names,param.names,drop=FALSE],center=VDATA$center[param.names]));
  return(mat.scaled[param.names,param.names,drop=FALSE]);
}

cfg.hzar.default.mcmc <- hzar.make.mcmcParam(chainLength=1e6,
                                         burnin=1e4,
                                         verbosity=5e5,
                                         thin=100);

cfg.hzar.quiet.mcmc   <- hzar.make.mcmcParam(chainLength=1e6,
                                         burnin=1e4,
                                         verbosity=0,
                                         thin=100);

##refitting
hzar.cov.mcmc<-function(clineLLfunc,mcmcRaw,pDiv=15,random=1e4,passCenter=FALSE){
  mcmc.nm<-colnames(mcmcRaw);
  pL<-lapply(mcmc.nm,function(x){min(mcmcRaw[,x])});
  names(pL)<-mcmc.nm;
  pU<-lapply(mcmc.nm,function(x){max(mcmcRaw[,x])});
  names(pU)<-mcmc.nm;
  return(hzar.cov.rect(clineLLfunc,pL,pU,pDiv,random,passCenter));
}
hzar.next.fitRequest <- function(oldFitRequest){
  seedChannel<-1;
  baseSeed=rep(12345,6)
  if(is.list(oldFitRequest$mcmcParam$seed)){
    seedChannel=oldFitRequest$mcmcParam$seed[[2]];
    baseSeed=oldFitRequest$mcmcParam$seed[[1]]
  }
  if(identical( attr(oldFitRequest,"fit.run") , TRUE)){
    seedChannel<-seedChannel+1;
  } else {
    seedChannel<-seedChannel+10;
  }
  mcmcParam=hzar.make.mcmcParam(oldFitRequest$mcmcParam$chainLength,
    oldFitRequest$mcmcParam$burnin,
    oldFitRequest$mcmcParam$verbosity,
    oldFitRequest$mcmcParam$thin,
    seedChannel,lecuyerSeed=baseSeed);
  mdlParam<-oldFitRequest$modelParam;
  covMatrix<-oldFitRequest$cM
  if(identical( attr(oldFitRequest,"fit.success") , TRUE)){
    nGen <- dim(oldFitRequest$mcmcRaw)[[1]]
    nSam <- min(nGen,5e3)
    mcmcSubset<-oldFitRequest$mcmcRaw[sample(nGen),];
    subLL<-hzar.eval.clineLL(mcmcSubset,oldFitRequest$llFunc);
    covData <- NULL
    try(covData <- cov.wt(mcmcSubset[subLL>max(subLL-4),]))
    if(sum(unique(subLL)>max(subLL-4))<0.95*nGen){
      try({
        ## if(sum(unique(subLL)>max(subLL-4))>1e3){
        ##   covData<-hzar.cov.mcmc(oldFitRequest$llFunc,mcmcSubset[subLL>max(subLL-4),],passCenter=TRUE);
        ## }else {
        ##   mcmcSubset<-oldFitRequest$mcmcRaw[sample(nGen),];
        ##   subLL<-hzar.eval.clineLL(mcmcSubset,oldFitRequest$llFunc);
        if(sum(unique(subLL)>max(subLL-4))>1e3){
          wt <- subLL[subLL>max(subLL-4),]
          wt <- exp(wt-max(wt))
          covData <- cov.wt(mcmcSubset[subLL>max(subLL-4),],wt)
        }else if(sum(unique(subLL)>max(subLL-4))>1e2){
          covData<-hzar.cov.mcmc(oldFitRequest$llFunc,mcmcSubset[subLL>max(subLL-4),],passCenter=TRUE);
        }else if(sum(unique(subLL)>max(subLL-10))>1e2){
          covData<-hzar.cov.mcmc(oldFitRequest$llFunc,mcmcSubset[subLL>max(subLL-10),],passCenter=TRUE);
        }else {
          covData<-hzar.cov.mcmc(oldFitRequest$llFunc,mcmcSubset,passCenter=TRUE);
        }
      },silent=TRUE)
    }
    if(!is.null(covData)){
      covMatrix<-covData$cov;
      new.center<-covData$center[names(mdlParam$init)];
      if(oldFitRequest$llFunc(new.center)>-1e6)
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
      try({ junk <-  solve(-naiveHessian(mdlParam$init,oldFitRequest$llFunc));
            junk <- appHScale(junk,mdlParam$lower, mdlParam$upper);
            if(all(diag(junk)>0)){
              chol(junk);
              covMatrix <- junk;}
          },silent=TRUE)
    }
    if(is.null(covMatrix)){
      try({
        junk <- appThetaWalkerR(mdlParam$init,
                                oldFitRequest$llFunc,
                                mdlParam$lower, 
                                mdlParam$upper, random = 1000,
                                passCenter=TRUE);
        covMatrix <- junk$cov
        mdlParam$init  <- junk$center
        
      })
    }
  }
  return(hzar.make.fitRequest(mdlParam,
                              covMatrix,
                              oldFitRequest$llFunc,
                              mcmcParam));
}



freqCompileLLF <- function(pObs,nObs,pExp){
  return(substitute((1-pObs)*log((1-pEst)/(1-pObs))*nObs+pObs*log(pEst/pObs)*nObs,
                    list(pObs=pObs,nObs=nObs,pEst=pExp)));
}

freqCompileLLEdge <- function(nObs,pExp){
  return(substitute(nObs*log(pEst),list(nObs=nObs,pEst=pExp)));
}
freqCompilePDF <- function(dist,pExp,target){
  bquote(.(target) <- .(eval(substitute(substitute(A,list(x=dist)),list(A=pExp)))))
}
#s.LIB=list(paren=
s.exp.LIB<-list()
s.exp.SYM <- c();
## sym=list(),use=list(),run=list()
s.exp.LIB.add <- function(sym,use,run){
  s.exp.LIB<<-c(s.exp.LIB,list(list(sym=c(sym),use=use,run=run)))
  s.exp.SYM<<-unique(c(s.exp.SYM,sym))
}

s.exp2NumArg <- function(x) length(x)==3 && is.numeric(x[[2]]) && is.numeric(x[[3]])
s.exp1NumArg <- function(x) length(x)==2 && is.numeric(x[[2]])
s.expNumSimpA <- function(x) length(x)==3 && is.numeric(x[[2]]) && length(x[[2]])==1
s.expNumSimpB <- function(x) length(x)==3 && is.numeric(x[[3]]) && length(x[[3]])==1
s.exp1LangArg <- function(x) length(x)==2 && !is.symbol(x[[2]]) && is.language(x[[2]])

s.exp2boolArg <- function(x) length(x)==3 && is.logical(x[[2]]) && is.logical(x[[3]])
s.exp1boolArg <- function(x) length(x)==2 && is.logical(x[[2]])
s.expBoolSimpA <- function(x) length(x)==3 && is.logical(x[[2]]) && length(x[[2]])==1
s.expBoolSimpB <- function(x) length(x)==3 && is.logical(x[[3]]) && length(x[[3]])==1

s.exp.LIB.add(sym=as.name("+"),use=s.exp2NumArg,run=function(x) x[[2]]+x[[3]])
s.exp.LIB.add(sym=as.name("-"),use=s.exp2NumArg,run=function(x) x[[2]]-x[[3]])
s.exp.LIB.add(sym=as.name("-"),use=s.exp1NumArg,run=function(x) -x[[2]])
s.exp.LIB.add(sym=as.name("-"),
              use=s.exp1LangArg,
              run=function(x) {
                if(x[[2]][[1]]==as.name("*")||x[[2]][[1]]==as.name("/")){
                  if(is.numeric(x[[2]][[2]])) {
                    x[[2]][[2]]=-x[[2]][[2]];
                    return(x[[2]])
                  }
                  if(is.numeric(x[[2]][[3]])) {
                    x[[2]][[3]]=-x[[2]][[3]];
                    return(x[[2]])
                  }
                }
                if(x[[2]][[1]]==as.name("-")){
                  if(length(x[[2]])==2)
                    return(x[[2]][[2]])
                  junk=x[[2]][[3]]
                  x[[2]][[3]]=x[[2]][[2]]
                  x[[2]][[2]]=junk
                  return(x[[2]])
                }
                if(x[[2]][[1]]==as.name("(")&&x[[2]][[2]][[1]]==as.name("-")){
                  if(length(x[[2]][[2]])==2)
                    return(x[[2]][[2]][[2]])
                  junk=x[[2]][[2]][[3]]
                  x[[2]][[2]][[3]]=x[[2]][[2]][[2]]
                  x[[2]][[2]][[2]]=junk
                  return(x[[2]][[2]])
                }
                x})
s.exp.LIB.add(sym=as.name("*"),use=s.exp2NumArg,run=function(x) x[[2]]*x[[3]])
s.exp.LIB.add(sym=as.name("/"),use=s.exp2NumArg,run=function(x) x[[2]]/x[[3]])

s.exp.LIB.add(sym=as.name("+"),use=s.expNumSimpA,run=function(x) {if(x[[2]]==0) return(x[[3]]); x})
s.exp.LIB.add(sym=as.name("-"),use=s.expNumSimpA,run=function(x) {if(x[[2]]==0) return(c(as.name("-"),x[[3]])); x})
s.exp.LIB.add(sym=c(as.name("+"),as.name("-")),use=s.expNumSimpB,run=function(x) {if(x[[3]]==0) return(x[[2]]); x})

s.exp.LIB.add(sym=as.name("*"),
              use=s.expNumSimpA,
              run=function(x) {
                if(x[[2]]==0) return(0);
                if(x[[2]]==1) return(x[[3]]);
                if(!is.symbol(x[[3]])&&is.language(x[[3]])){
                  if((x[[3]][[1]]==as.name("*")||x[[3]][[1]]==as.name("/"))
                     &&is.numeric(x[[3]][[2]])) {
                    x[[3]][[2]]=x[[2]]*x[[3]][[2]];
                    return(x[[3]])
                  }
                  if(x[[3]][[1]]==as.name("*")&&is.numeric(x[[3]][[3]])) {
                    x[[3]][[3]]=x[[2]]*x[[3]][[3]];
                    return(x[[3]])
                  }
                }
                x})
## s.exp.LIB.add(sym=as.name("*"),use=s.expNumSimpA,run=function(x) {if(x[[2]]==1) return(x[[3]]); x})
s.exp.LIB.add(sym=as.name("*"),use=s.expNumSimpB,run=function(x) {if(x[[3]]==0) return(0); x})
s.exp.LIB.add(sym=c(as.name("*"),as.name("/")),
              use=s.expNumSimpB,
              run=function(x) {
                if(x[[3]]==1) return(x[[2]]);
                if(x[[3]]==-1) return(c(as.name("-"),x[[2]]));
                x})

s.exp.LIB.add(sym=as.name("!"),use=s.exp1boolArg,run=function(x) !x[[2]])
s.exp.LIB.add(sym=as.name("|"),use=s.exp2boolArg,run=function(x) x[[2]]|x[[3]])
s.exp.LIB.add(sym=as.name("&"),use=s.exp2boolArg,run=function(x) x[[2]]&x[[3]])
s.exp.LIB.add(sym=c(as.name("&&"),as.name("&")),use=s.expBoolSimpA,run=function(x) {if(x[[2]])return(x[[3]]);FALSE})
s.exp.LIB.add(sym=c(as.name("||"),as.name("|")),use=s.expBoolSimpA,run=function(x) {if(!x[[2]])return(x[[3]]);TRUE})
s.exp.LIB.add(sym=c(as.name("&&"),as.name("&")),use=s.expBoolSimpB,run=function(x) {if(x[[3]])return(x[[2]]);FALSE})
s.exp.LIB.add(sym=c(as.name("||"),as.name("|")),use=s.expBoolSimpB,run=function(x) {if(!x[[3]])return(x[[2]]);TRUE})
s.exp.LIB.add(sym=as.name("<"),use=s.exp2NumArg,run=function(x) x[[2]]<x[[3]])
s.exp.LIB.add(sym=as.name(">"),use=s.exp2NumArg,run=function(x) x[[2]]>x[[3]])
s.exp.LIB.add(sym=as.name(">="),use=s.exp2NumArg,run=function(x) x[[2]]>=x[[3]])
s.exp.LIB.add(sym=as.name("<="),use=s.exp2NumArg,run=function(x) x[[2]]<=x[[3]])
s.exp.LIB.add(sym=as.name("=="),use=s.exp2NumArg,run=function(x) x[[2]]==x[[3]])
s.exp.LIB.add(sym=as.name("!="),use=s.exp2NumArg,run=function(x) x[[2]]!=x[[3]])


s.exp.LIB.add(
  sym=as.name("("),
  use=function(x) length(x)==2 && x[[1]]==as.name("(") && (is.symbol(x[[2]]) || !is.language(x[[2]])),
  run=function(x) x[[2]])


simplify.exp <- function(expL){
  ##cat("l")
  if(is.symbol(expL)||is.numeric(expL)||is.logical(expL))return(expL)
  ##cat("L")
  proposed <- lapply(expL,simplify.exp)
  ##cat("\ns")
  if(is.symbol(proposed[[1]])&&all(proposed[1]%in% s.exp.SYM )){
    for(iter in 1:length(s.exp.LIB)){
      
      ##cat("_")
      layer=s.exp.LIB[[iter]]
      
      if(is.symbol(proposed[[1]])&&all(proposed[1]%in% layer$sym) && layer$use(proposed)){
        ##cat("+")
        ##str(proposed)
        proposed <- layer$run(proposed)
      }
      if(is.symbol(proposed))return(proposed)
      if(is.language(proposed))proposed <- as.list(proposed)
      if(!is.list(proposed))return(proposed)
    }
  }
  return(as.call(proposed))
}

#s.LIB=list(paren=
rA.exp.LIB<-list()
rA.exp.SYM <- c();
## sym=list(),use=list(),run=list()
rA.exp.LIB.add <- function(sym,use,run){
  rA.exp.LIB<<-c(rA.exp.LIB,list(list(sym=c(sym),use=use,run=run)))
  rA.exp.SYM<<-unique(c(rA.exp.SYM,sym))
}

rA.exp2Arg <- function(x) length(x)==3
rA.exp.LIB.add(sym=c(as.name(">"),as.name(">=")),
               use=rA.exp2Arg,
               run=function(x) bquote(ifelse(
                 .(x[[2]]) %lt% .(x[[3]]),
                 .(x[[2]]) -.(x[[3]])-1e6,
                 0)))
rA.exp.LIB.add(sym=c(as.name("<"),as.name("<=")),
               use=rA.exp2Arg,
               run=function(x) bquote(ifelse(
                 .(x[[2]]) %gt% .(x[[3]]),
                 .(x[[3]])-1e6 -.(x[[2]]),
                 0)))
rA.exp.LIB.add(sym=as.name("|"),use=rA.exp2Arg,run=function(x) bquote (.(x[[2]])*.(x[[3]])))
rA.exp.LIB.add(sym=as.name("&"),use=rA.exp2Arg,run=function(x) bquote (.(x[[2]])+.(x[[3]])))

#s.LIB=list(paren=
rB.exp.LIB<-list()
rB.exp.SYM <- c();
## sym=list(),use=list(),run=list()
rB.exp.LIB.add <- function(sym,use,run){
  rB.exp.LIB<<-c(rB.exp.LIB,list(list(sym=c(sym),use=use,run=run)))
  rB.exp.SYM<<-unique(c(rB.exp.SYM,sym))
}

rB.exp2Arg <- function(x) length(x)==3
rB.exp.LIB.add(sym=as.name("%gt%"),
               use=rB.exp2Arg,
               run=function(x) bquote(.(x[[2]]) > .(x[[3]]) ))

rB.exp.LIB.add(sym=as.name("%lt%"),
               use=rB.exp2Arg,
               run=function(x) bquote(.(x[[2]]) < .(x[[3]]) ))
               

req.Expand.exp <- function(expL){
  if(is.symbol(expL)||is.numeric(expL)||is.logical(expL))return(expL)
  proposed <- lapply(expL,req.Expand.exp)
  if(is.symbol(proposed[[1]])&&all(proposed[1]%in% rA.exp.SYM )){
    for(iter in 1:length(rA.exp.LIB)){
      layer=rA.exp.LIB[[iter]]
      
      if(is.symbol(proposed[[1]])&&all(proposed[1]%in% layer$sym) && layer$use(proposed)){
        ##cat("+")
        ##str(proposed)
        proposed <- layer$run(proposed)
      }
      if(is.symbol(proposed))return(proposed)
      if(is.language(proposed))proposed <- as.list(proposed)
      if(!is.list(proposed))return(proposed)
    }
  }
  return(as.call(proposed))
}

req.Collapse.exp <- function(expL){
  if(is.symbol(expL)||is.numeric(expL)||is.logical(expL))return(expL)
  proposed <- lapply(expL,req.Collapse.exp)
  if(is.symbol(proposed[[1]])&&all(proposed[1]%in% rB.exp.SYM )){
    for(iter in 1:length(rB.exp.LIB)){
      layer=rB.exp.LIB[[iter]]
      
      if(is.symbol(proposed[[1]])&&all(proposed[1]%in% layer$sym) && layer$use(proposed)){
        ##cat("+")
        ##str(proposed)
        proposed <- layer$run(proposed)
      }
      if(is.symbol(proposed))return(proposed)
      if(is.language(proposed))proposed <- as.list(proposed)
      if(!is.list(proposed))return(proposed)
    }
  }
  return(as.call(proposed))
}

req.edge.exp <- function(expL) simplify.exp(req.Collapse.exp(req.Expand.exp(expL)))



freq.LLfunc <- function(obsData, model,tInit,tFixed,
                        LLrejectedModel = -1e+08){
    baseFunc <- function(theta) 0;
    model.req=model$req;
    model.gen=model$func
    model.obsData=obsData
    param.fixed=tFixed
    frame=obsData$frame
    frameA=subset(frame,frame$obsFreq == 0)
    frameB=subset(frame,frame$obsFreq == 1)
    pExp=model$pExp
    pExp <- ll.compile.theta(tInit,tFixed,pExp)
    ## print(pExp)
    gLLc<-list()
    pDF <- list()
    frame=subset(frame,frame$obsFreq !=0 & frame$obsFreq != 1)
    new.formals=c(formals(model.req)[names(tInit)],tFixed)
    ## new.formals=tInit
    if(nrow(frame) > 1){
      pDF <- c(pDF,freqCompilePDF(quote(frame$dist), pExp,quote(pFunc)))
      gLLc <- c(gLLc,freqCompileLLF(quote(frame$obsFreq),quote(frame$n), quote(pFunc)))
    }else if(nrow(frame) == 1){
      pDF <- c(pDF,freqCompilePDF(frame$dist, pExp,quote(pFunc)))
      gLLc <- c(gLLc,freqCompileLLF(frame$obsFreq,frame$n, quote(pFunc)))
    }
    if(nrow(frameA)>1){
      qExp=bquote(1-.(pExp))
      pDF <- c(pDF,freqCompilePDF(quote(frameA$dist), qExp,quote(qEdge)))
      gLLc <- c(gLLc,freqCompileLLEdge(quote(frameA$n), quote(qEdge)))
    }else if(nrow(frameA)==1){
      qExp=bquote(1-.(pExp))
      pDF <- c(pDF,freqCompilePDF(frameA$dist, qExp,quote(qEdge)))
      gLLc <- c(gLLc,freqCompileLLEdge(frameA$n, quote(qEdge)))
    }
    if(nrow(frameB)>1){
      pDF <- c(pDF,freqCompilePDF(quote(frameB$dist), pExp,quote(pEdge)))
      gLLc <- c(gLLc,freqCompileLLEdge(quote(frameB$n), quote(pEdge)))
    }else if(nrow(frameB)==1){
      pDF <- c(pDF,freqCompilePDF(frameB$dist, pExp,quote(pEdge)))
      gLLc <- c(gLLc,freqCompileLLEdge(frameB$n, quote(pEdge)))
    }
    ## print(pDF)
    ## print(gLLc)
    if(length(gLLc)==0) stop("No observed data?")
    gLL <- as.call(c(quote(sum),gLLc))
    ## print(gLL)
    reqL <- ll.compile.theta(tInit,tFixed,body(model$req)[[2]])
    LLrejectedModel <- bquote(.(LLrejectedModel)+.(req.edge.exp(reqL)))
    gLL <-  step1VectorExpF(reqL,
                            gLL,
                            LLrejectedModel)
    ## print(gLL)
    ## gLL <- eval(substitute(substitute(LLfunc,tF),list(LLfunc=gLL,tF=tFixed)))
    ##print(gLL)
    body(baseFunc) <-simplify.exp(as.call(c(as.name("{"),
                               pDF,
                               bquote(res <- .(gLL)),
                               expression(if(any(is.na(res))) print(theta)),
                               bquote(ifelse(is.na(res),.(LLrejectedModel),res)))));
    llFunc=baseFunc
    old.formals=formals(model.req)
    formals(model.req) <- c(old.formals[names(tInit)],tFixed)
    formals(model.gen) <- c(old.formals[names(tInit)],tFixed)
    ## eval(substitute(substitute(LLfunc,eL),
    ##                                 list(LLfunc=gLL,eL=tMap)))
    ## body(baseFunc) <- substitute(evalq(LLfunc,envir=theta),list(LLfunc=gLL))
    ## environment(baseFunc) <- .GlobalEnv 
    baseFunc
}
hzar.first.fitRequest.old.ML <-function(model,obsData,verbose=TRUE){
  
  if(verbose){
    mcmcParam<-cfg.hzar.default.mcmc;
  }else {
    mcmcParam<-cfg.hzar.quiet.mcmc;
  }
   ## print("A");
  
  modelParam<-splitParameters(model$parameterTypes);
   ## print("A");
  clineLLfunc <- freq.LLfunc(obsData,model,modelParam$init,modelParam$fixed)
  ## clineLLfunc<-hzar.make.clineLLfunc.old.ML(names(modelParam$init),
##                                             modelParam$fixed,
##                                             model$req,
##                                             model$func,
##                                             obsData$model.LL);
   ## print("A");
  covMatrix<-NULL;
  try(  covMatrix<-hzar.cov.rect(clineLLfunc,modelParam$lower,modelParam$upper,random=1e4));
   ## print("A");
  return(hzar.make.fitRequest(modelParam,covMatrix,clineLLfunc,mcmcParam));
}


## hzar.multiFit.doNext

hzar.chain.doSeq <- function(hzar.request, count=3,collapse=FALSE,announce.complete="Chain Complete"){
  if(collapse){
    mcmcParam=hzar.request$mcmcParam;
    mcmcParam$chainLength= mcmcParam$chainLength*count;
    mcmcParam$burnin= mcmcParam$burnin*count;
    mdlParam=hzar.request$modelParam;
    
  }
  hzar.results<-list();
  for(iter in 1:count){
    print(iter);
    if(iter>1){
      hzar.request<-NULL;
      print(summary(try(hzar.request<-
                        hzar.next.fitRequest(hzar.results[[iter-1]] ))));
    }
    if(!inherits(hzar.request,c("hzar.fitRequest"))){
      warning("Failed to generate next fit request. Returning successful runs.");
      
      return(hzar.results);
    }
    hzar.results[[iter]] <- hzar.request;
    ##format(hzar.request$cM,width=8);
    print(summary(try(hzar.results[[iter]] <- hzar.doFit(hzar.request))));
    ##print(attributes(hzar.results[[iter]]))
  }
  print(announce.complete);
  if(collapse){
    raw.data<-do.call(rbind,lapply(hzar.results,function(x)x$mcmcRaw));
    
    rawMCMC<-mcmc(data=raw.data,thin=thin(hzar.results[[count]]$mcmcRaw),start=mcmcParam$burnin+1);
    names(rawMCMC)<-names(hzar.results[[count]]$mcmcRaw);
    return(list(hzar.make.fitRequest(mdlParam,
                                     hzar.results[[count]]$cM,
                                     hzar.results[[count]]$llFunc,
                                     mcmcParam,
                                     mcmcRaw=rawMCMC,
                                     TRUE,
                                     prod(as.logical(lapply(hzar.results,attr,"fit.success"))))));
  }
  return(hzar.results);
}
                                                                                                              
