## So I need this file to handle the setup of the data and the model.

## First, i will define the sample level likelihood Function.

sampleLikelihoodMolecularPop=function(pEst,pObs,N) {
  p1=N*log(pEst);
  p0=N*log(1-pEst);
  lP=((1-pObs)*log((1-pEst)/(1-pObs))+pObs*log(pEst/pObs))*N;
  tmp=ifelse(pObs==1,p1,lP);
  res=ifelse(pObs==0,p0,tmp);
  return(res);
}

## Second, i will setup the data.
hzar.doMolecularData1DPops<-function(distance,pObs,nEff,
                                     siteID=paste("P",1:length(distance),sep=""),
                                     ylim=extendrange(c(0,1))){
  if((length(distance) != length(pObs))     ||
     (length(distance) != length(nEff)) ||
     (length(pObs) != length(nEff))  ){
    stop("Distance, pObs and nEff are not of the same length.");
  }
  obj<-list(frame=na.omit(data.frame(
              dist=distance,
              obsFreq=pObs,
              n=nEff,
              row.names=siteID)));
  obj$ylim <- ylim;
  obj$model.LL <- function(model.func){
    pEst=model.func(obj$frame$dist);
##res<-numeric(length(pEst));
  ##  for(iter in 1:(length(pEst))){
      res<-sampleLikelihoodMolecularPop(pEst=as.numeric(pEst),
                                        pObs=as.numeric(obj$frame$obsFreq),
                                        N=as.numeric(obj$frame$n));
    ##}
                                                
    result<-sum(res);
    if(is.na(result))
      return(-1e8);
    return(result);
  }
class(obj)<-c("clineSampleData1D","hzar.obsData");
  return(obj);
}

## TODO: Add additional parameter defaults, including d1
## (center-deltaL), d2 (deltaR-center), betaL and betaR (rather
## complex).

##format of parameter data frame?
mkParam<-function(name,value,weight,lower,upper,isFixed=FALSE,isRB01=FALSE){
  objCParamC<-list(val=value,w=weight);
  class(objCParamC)<-"clineParameter";
  attr(objCParamC,"param")<-name
  attr(objCParamC,"fixed")<-isFixed;
  attr(objCParamC,"limit.lower")<- lower;
  attr(objCParamC,"limit.upper")<- upper;
  attr(objCParamC,"realBTWN01")<- isRB01
  return(objCParamC);
}
CLINEPARAMETERS<-list(center=mkParam("center",10,1.5,-1e8,1e8),
                      width =mkParam("width", 10,1.5,1e-8,1e8),
                      pMin  =mkParam("pMin",  0, 1.1,   0,  1,isRB01=TRUE),
                      pMax  =mkParam("pMax",  1, 1.1,   0,  1,isRB01=TRUE),
                      xMin  =mkParam("xMin",0,1.5,-1e8,1e8),
                      xMax  =mkParam("xMax",10,1.5,-1e8,1e8),
                      muL=mkParam("muL",1,1.5,  -1e8, 1e8),
                      muR=mkParam("muR",1,1.5,  -1e8, 1e8),
                      varL=mkParam("varL",1,1.5,  1e-8, 1e8),
                      varR=mkParam("varR",1,1.5,  1e-8, 1e8),
                      varH=mkParam("varH",1,1.5,  1e-8, 1e8),
                      deltaL=mkParam("deltaL",1,1.5,  1e-8, 1e8),
                      deltaR=mkParam("deltaR",1,1.5,  1e-8, 1e8),
                      deltaM=mkParam("deltaM",1,1.5,  1e-8, 1e8),
                      tauL  =mkParam("tauL", 0.5,1.1,  0,   1,isRB01=TRUE),
                      tauR  =mkParam("tauR", 0.5,1.1,  0,   1,isRB01=TRUE),
                      tauM  =mkParam("tauM", 0.5,1.1,  0,   1,isRB01=TRUE));

## Internal getter and setter methods for parameters 

param.check <- function(x){
  if(!inherits(x,"clineParameter"))
    stop("Argument x is not a clineParameter")
}
param.init <- function(x){ param.check(x); x$val }
"param.init<-" <- function(x,value){ param.check(x); x$val<-value; x }
param.tune <- function(x){ param.check(x); x$w }
"param.tune<-" <- function(x,value){ param.check(x); x$w<-value; x }
param.fix <- function(x){ param.check(x); attr(x,"fixed") }
"param.fix<-" <- function(x,value){ param.check(x);  attr(x,"fixed")<-as.logical(value); x }
param.upper <- function(x){ param.check(x); attr(x,"limit.upper") }
"param.upper<-" <- function(x,value){ param.check(x);  attr(x,"limit.upper")<-value; x }
param.lower <- function(x){ param.check(x); attr(x,"limit.lower") }
"param.lower<-" <- function(x,value){ param.check(x);  attr(x,"limit.lower")<-value; x }
param.name <- function(x){ param.check(x); attr(x,"param") }
print.clineParameter <- function(x,...){
  print(data.frame(param=param.name(x),
                   init=param.init(x),
                   tune=param.tune(x),
                   fixed=param.fix(x),
                   lower=param.lower(x),
                   upper=param.upper(x)));
  invisible(x);
}

## suggest upper and lower bounds for cov matrix based on observed data:
## have vector of distances (obsData$frame$dist),
## and a vector of data points (obsData$frame$pObs?).

cline.suggestionFunc1D <-
  list(center = function(dist,obs){
    return(c(min(dist),max(dist)));},
       width  = function(dist,obs){
         return(c(0,max(dist)-min(dist)));},
       pMin   = function(dist,obs){return(c(0,1));},
       pMax   = function(dist,obs){return(c(0,1));},
       xMin   = function(dist,obs){
         delta=(max(obs)-min(obs));
         return(c(min(obs)-delta/10,max(obs)+delta/10));},
       xMax   = function(dist,obs){
         delta=(max(obs)-min(obs));
         return(c(min(obs)-delta/10,max(obs)+delta/10));},
       deltaL = function(dist,obs){
         return(c(0,max(dist)-min(dist)));},
       deltaR = function(dist,obs){
         return(c(0,max(dist)-min(dist)));},
       deltaM = function(dist,obs){
         return(c(0,max(dist)-min(dist)));},
       tauL   = function(dist,obs){return(c(0,1));},
       tauR   = function(dist,obs){return(c(0,1));},
       tauM   = function(dist,obs){return(c(0,1));});
       
      
## I need a revised cline meta model, that is less memory instensive
## and creates less function objects.  It also needs sufficient
## modularity that imposing arbitrary priors, using a different
## parameterization, or compositing models of multiple loci (or
## traits) is feasible.


## ## format of the cline meta model
## objCMeta<-list();
## objCMeta$parameterTypes<-list("center","width");
## objCMeta$req<-function(center,width);
## objCMeta$func<-function(center,width);
## class(objCMeta)<-"clineMetaModel";

## Internal dispatch methods for handling parameters


meta.check <- function(x){
  if(!inherits(x,"clineMetaModel"))
    stop("Argument x is not a clineMetaModel")
}

meta.param.names <- function(x) names(x$parameterTypes)
meta.check.names <- function(x,cVec) {
  meta.param.names(x)-> pNames;
  if(length(pNames)!=length(cVec)
     || !all(pNames %in% cVec)
     || !all(cVec %in% pNames) )
    stop(paste("Invalid Assignment of names:", paste(cVec,collapse=", "),"\n"))
  return(TRUE)
}

meta.fixValue <- function(x,value,tC=is.numeric){
  if(tC(value)&&length(value)==1){
    value<-as.list(rep(value,length(x$parameterTypes)))
    names(value) <- meta.param.names(x);
  }else{
    if(!is.list(value)) value <- as.list(value)
    meta.check.names(x,names(value))
  }
  value
}
meta.param.m <- function(x,pF){
  res <- lapply(x$parameterTypes,pF);
  names(res)<-meta.param.names(x);
  res
}

 

meta.init <- function(x) meta.param.m(x,param.init)
"meta.init<-" <- function(x,value) {
  value <- meta.fixValue(x,value)
  for(iter in meta.param.names(x))
    param.init(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 

meta.tune <- function(x) meta.param.m(x,param.tune)
"meta.tune<-" <- function(x,value) {
  value <- meta.fixValue(x,value)
  for(iter in meta.param.names(x))
    param.tune(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 

meta.fix <- function(x) meta.param.m(x,param.fix)
"meta.fix<-" <- function(x,value) { 
  value <- meta.fixValue(x,value,is.logical)
  for(iter in meta.param.names(x))
    param.fix(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 

meta.lower <- function(x) meta.param.m(x,param.lower)
"meta.lower<-" <- function(x,value) { 
  value <- meta.fixValue(x,value)
  for(iter in meta.param.names(x))
    param.lower(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 

meta.upper <- function(x) meta.param.m(x,param.upper)
"meta.upper<-" <- function(x,value) {
  value <- meta.fixValue(x,value)
  for(iter in meta.param.names(x))
    param.upper(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 



hzar.meta.init <- function(x) meta.param.m(x,param.init)
"hzar.meta.init<-" <- function(x,value) {
  value <- meta.fixValue(x,value)
  for(iter in meta.param.names(x))
    param.init(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 

hzar.meta.tune <- function(x) meta.param.m(x,param.tune)
"hzar.meta.tune<-" <- function(x,value) {
  value <- meta.fixValue(x,value)
  for(iter in meta.param.names(x))
    param.tune(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 

hzar.meta.fix <- function(x) meta.param.m(x,param.fix)
"hzar.meta.fix<-" <- function(x,value) { 
  value <- meta.fixValue(x,value,is.logical)
  for(iter in meta.param.names(x))
    param.fix(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 

hzar.meta.lower <- function(x) meta.param.m(x,param.lower)
"hzar.meta.lower<-" <- function(x,value) { 
  value <- meta.fixValue(x,value)
  for(iter in meta.param.names(x))
    param.lower(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 

hzar.meta.upper <- function(x) meta.param.m(x,param.upper)
"hzar.meta.upper<-" <- function(x,value) {
  value <- meta.fixValue(x,value)
  for(iter in meta.param.names(x))
    param.upper(x$parameterTypes[[iter]]) <- value[[iter]];
  x
} 


print.clineMetaModel <- function(x,...){
  for(iter in names(x)){
    if(iter!="parameterTypes"){
      cat(paste(iter,":\n",sep=""));
      print(x[[iter]],...)
    }
  }
  cat("Cline Parameters:\n")
  print(data.frame(row.names=meta.param.names(x),
                   init=as.numeric(meta.init(x)),
                   tune=as.numeric(meta.tune(x)),
                   fixed=as.logical(meta.fix(x)),
                   lower=as.numeric(meta.lower(x)),
                   upper=as.numeric(meta.upper(x))),...);
  invisible(x);
}



cline.func.gen <- function(model){
  f <- model$req
  pExp <- model$pExp
  bF <- bquote(return(function(x) .(pExp)))
  body(f) <- bF
##  attr(f,"source") <- deparse(bF)
  environment(f) <- .GlobalEnv
  model$func <- f
  model
}
cline.u.ascend <- quote((x - center) * 4/width)
cline.du.ascend <- quote(x - center)
cline.u.descend <- quote((x - center) * -4/width)
cline.du.descend <- quote(center-x)
cline.meta.simple.scaled.ascending =
  cline.func.gen(
  list(
       prior=function(center,width,pMin,pMax){
  return(0); },
       pExp=cline.exp.scale(cline.exp.sigmoid(quote((x - center) * 4/width))),
       
      ## func=function(center,width,pMin,pMax){
##          pCline<- function(x) {
##            u <- (x - center) * 4/width
##            return(pMin+(pMax-pMin)* (1/(1+ exp(-u)))) }
         
##          return(pCline)
##        },
       req=function(center,width,pMin,pMax)
         return(pMin>=0 & pMax<=1 & pMin<pMax & width>0),
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax")]
       ));
cline.meta.simple.scaled.descending =
  cline.func.gen(list(
       prior=function(center,width,pMin,pMax){
  return(0); },
      pExp=cline.exp.scale(cline.exp.sigmoid(quote((x - center) * -4/width))),
       ## func=function(center,width,pMin,pMax){
##          pCline<- function(x) {
##            u <- (x - center) * -4/width
##            return(pMin+(pMax-pMin)* (1/(1+ exp(-u)))) }
         
##          return(pCline)
##        },
       req=function(center,width,pMin,pMax)
         return(pMin>=0 & pMax<=1 & pMin<pMax & width>0),
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax")]
       ));
class(cline.meta.simple.scaled.ascending)<-"clineMetaModel";
class(cline.meta.simple.scaled.descending)<-"clineMetaModel";
cline.meta.tailed.scaled.ascending =
  cline.func.gen(list(req= function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR)
       
         return(width>0 & deltaL>=0 & deltaR>=0 &
               (pMax-pMin)* (deltaL+deltaR)<=width*50 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauL>=0 & tauL<=1 &
                tauR>=0 & tauR<=1 ),
       prior=function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR){
         return(0);
         ## return(-2*log(width/2)-2*(deltaL+deltaR)/width);
       },
       pExp=cline.exp.scale(cline.exp.stepBoth(
         cline.u.ascend,cline.du.ascend,
         quote(deltaL),
         quote(deltaR),
         cline.exp.lower(cline.u.ascend,quote(deltaL),quote(tauL)),
         cline.exp.upper(cline.u.ascend,quote(deltaR),quote(tauR)))),
      ## func=function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR)
##        {
##          gamma=4/width;
##          tail.LO=meta.tail.lower(gamma=gamma,d1=deltaL,tau1=tauL);
##          tail.HI=meta.tail.upper(gamma=gamma,d2=deltaR,tau2=tauR);
##          clineComposite=
##            meta.cline.func.stepBoth(center=center,
##                                     direction=1,
##                                     gamma=gamma,
##                                     lowerTail=tail.LO,
##                                     upperTail=tail.HI);
##          return(meta.cline.func.pScale(pMin,pMax,clineComposite));
##        },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaL","tauL","deltaR","tauR")]
       ));
cline.meta.tailed.scaled.descending =
  cline.func.gen(list(req= function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR)
         return(width>0 & deltaL>=0 & deltaR>=0 &
                (pMax-pMin)*(deltaL+deltaR)<=width*50 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauL>=0 & tauL<=1 &
                tauR>=0 & tauR<=1 ),
       prior=function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR){
         return(0);
         ## return(-2*log(width/2)-2*(deltaL+deltaR)/width);
       },
       pExp=cline.exp.scale(cline.exp.stepBoth(
         cline.u.descend,cline.du.descend,
         quote(deltaR),
         quote(deltaL),
         cline.exp.lower(cline.u.descend,quote(deltaR),quote(tauR)),
         cline.exp.upper(cline.u.descend,quote(deltaL),quote(tauL)))),
      ## func=function(center,width,pMin,pMax,deltaL,tauL,deltaR,tauR)
##        {
##          gamma=4/width;
##          tail.LO=meta.tail.lower(gamma=gamma,d1=deltaR,tau1=tauR);
##          tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
##          clineComposite=
##            meta.cline.func.stepBoth(center=center,
##                                     direction=-1,
##                                     gamma=gamma,
##                                     lowerTail=tail.LO,
##                                     upperTail=tail.HI);
##          return(meta.cline.func.pScale(pMin,pMax,clineComposite));
##        },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaL","tauL","deltaR","tauR")]
       ));
class(cline.meta.tailed.scaled.ascending)<-"clineMetaModel";
class(cline.meta.tailed.scaled.descending)<-"clineMetaModel";

cline.meta.mtail.scaled.descending =
  cline.func.gen(list(req= function(center,width,pMin,pMax,deltaM,tauM)
         return(width>0 & deltaM>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauM>=0 & tauM<=1 ),
       prior=function(center,width,pMin,pMax,deltaM,tauM){
  return(0); },
       pExp=cline.exp.scale(cline.exp.stepBoth(
         cline.u.descend,cline.du.descend,
         quote(deltaM),
         quote(deltaM),
         cline.exp.lower(cline.u.descend,quote(deltaM),quote(tauM)),
         cline.exp.upper(cline.u.descend,quote(deltaM),quote(tauM)))),
      ## func=function(center,width,pMin,pMax,deltaM,tauM)
##        {
##          gamma=4/width;
##          tail.LO=meta.tail.lower(gamma=gamma,d1=deltaM,tau1=tauM);
##          tail.HI=meta.tail.upper(gamma=gamma,d2=deltaM,tau2=tauM);
##          clineComposite=
##            meta.cline.func.stepBoth(center=center,
##                                     direction=-1,
##                                     gamma=gamma,
##                                    lowerTail=tail.LO,
##                                     upperTail=tail.HI);
##          return(meta.cline.func.pScale(pMin,pMax,clineComposite));
##        },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaM","tauM")]
       ));
cline.meta.mtail.scaled.ascending =
  cline.func.gen(list(req= function(center,width,pMin,pMax,deltaM,tauM)
         return(width>0 & deltaM>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauM>=0 & tauM<=1 ),
       prior=function(center,width,pMin,pMax,deltaM,tauM){
  return(0); },
       pExp=cline.exp.scale(cline.exp.stepBoth(
         cline.u.ascend,cline.du.ascend,
         quote(deltaM),
         quote(deltaM),
         cline.exp.lower(cline.u.ascend,quote(deltaM),quote(tauM)),
         cline.exp.upper(cline.u.ascend,quote(deltaM),quote(tauM)))),
##       func=function(center,width,pMin,pMax,deltaM,tauM)
##        {
##          gamma=4/width;
##          tail.LO=meta.tail.lower(gamma=gamma,d1=deltaM,tau1=tauM);
##         tail.HI=meta.tail.upper(gamma=gamma,d2=deltaM,tau2=tauM);
##          clineComposite=
##            meta.cline.func.stepBoth(center=center,
##                                     direction=1,
##                                     gamma=gamma,
##                                     lowerTail=tail.LO,
##                                     upperTail=tail.HI);
##          return(meta.cline.func.pScale(pMin,pMax,clineComposite));
##        },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaM","tauM")]
       ));
class(cline.meta.mtail.scaled.ascending)<-"clineMetaModel";
class(cline.meta.mtail.scaled.descending)<-"clineMetaModel";
cline.meta.ltail.scaled.descending =
  cline.func.gen(list(req= function(center,width,pMin,pMax,deltaL,tauL)
         return(width>0 & deltaL>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauL>=0 & tauL<=1 ),
       prior=function(center,width,pMin,pMax,deltaL,tauL){
  return(0); },
       pExp=cline.exp.scale(cline.exp.stepUp(
         cline.u.descend,cline.du.descend,
         quote(deltaL),
         cline.exp.upper(cline.u.descend,quote(deltaL),quote(tauL)))),
      ## func=function(center,width,pMin,pMax,deltaL,tauL)
##        {
##          gamma=4/width;
##         # tail.LO=meta.tail.lower(gamma=gamma,d1=deltaL,tau1=tauL);
##          tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
##          clineComposite=
##            meta.cline.func.upStep(center=center,
##                                     direction=-1,
##                                     gamma=gamma,
##                               #      lowerTail=tail.LO,
##                                     upperTail=tail.HI);
##          return(meta.cline.func.pScale(pMin,pMax,clineComposite));
##        },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaL","tauL")]
       ));
cline.meta.ltail.scaled.ascending =
  cline.func.gen(list(req= function(center,width,pMin,pMax,deltaL,tauL)
         return(width>0 & deltaL>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauL>=0 & tauL<=1 ),
       prior=function(center,width,pMin,pMax,deltaL,tauL){
  return(0); },
       pExp=cline.exp.scale(cline.exp.stepLow(
         cline.u.ascend,cline.du.ascend,
         quote(deltaL),
         cline.exp.lower(cline.u.ascend,quote(deltaL),quote(tauL)))),
##       func=function(center,width,pMin,pMax,deltaL,tauL)
##        {
##          gamma=4/width;
##          tail.LO=meta.tail.lower(gamma=gamma,d1=deltaL,tau1=tauL);
##         # tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
##          clineComposite=
##            meta.cline.func.lowStep(center=center,
##                                     direction=1,
##                                     gamma=gamma,
##                                     lowerTail=tail.LO);
##          return(meta.cline.func.pScale(pMin,pMax,clineComposite));
##        },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaL","tauL")]
       ));
class(cline.meta.ltail.scaled.ascending)<-"clineMetaModel";
class(cline.meta.ltail.scaled.descending)<-"clineMetaModel";
cline.meta.rtail.scaled.ascending =
  cline.func.gen(list(req= function(center,width,pMin,pMax,deltaR,tauR)
         return(width>0 & deltaR>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauR>=0 & tauR<=1 ),
       prior=function(center,width,pMin,pMax,deltaR,tauR){
  return(0); },
       pExp=cline.exp.scale(cline.exp.stepUp(
         cline.u.ascend,cline.du.ascend,
         quote(deltaR),
         cline.exp.upper(cline.u.ascend,quote(deltaR),quote(tauR)))),
##       func=function(center,width,pMin,pMax,deltaR,tauR)
##        {
##          gamma=4/width;
##         # tail.LO=meta.tail.lower(gamma=gamma,d1=deltaL,tau1=tauL);
##          tail.HI=meta.tail.upper(gamma=gamma,d2=deltaR,tau2=tauR);
##          clineComposite=
##            meta.cline.func.upStep(center=center,
##                                     direction=1,
##                                     gamma=gamma,
##                               #      lowerTail=tail.LO,
##                                     upperTail=tail.HI);
##          return(meta.cline.func.pScale(pMin,pMax,clineComposite));
##        },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaR","tauR")]
       ));
cline.meta.rtail.scaled.descending =
  cline.func.gen(list(req= function(center,width,pMin,pMax,deltaR,tauR)
         return(width>0 & deltaR>=0 &
                pMin>=0 & pMax<=1 & pMin<pMax & 
                tauR>=0 & tauR<=1 ),
       prior=function(center,width,pMin,pMax,deltaR,tauR){
  return(0); },
       pExp=cline.exp.scale(cline.exp.stepLow(
         cline.u.descend,cline.du.descend,
         quote(deltaR),
         cline.exp.lower(cline.u.descend,quote(deltaR),quote(tauR)))),
##       func=function(center,width,pMin,pMax,deltaR,tauR)
##        {
##          gamma=4/width;
##          tail.LO=meta.tail.lower(gamma=gamma,d1=deltaR,tau1=tauR);
##         # tail.HI=meta.tail.upper(gamma=gamma,d2=deltaL,tau2=tauL);
##          clineComposite=
##            meta.cline.func.lowStep(center=center,
##                                     direction=-1,
##                                     gamma=gamma,
##                                     lowerTail=tail.LO);
##          return(meta.cline.func.pScale(pMin,pMax,clineComposite));
##        },
       parameterTypes=CLINEPARAMETERS[c("center","width","pMin","pMax","deltaR","tauR")]
       ));
class(cline.meta.rtail.scaled.ascending)<-"clineMetaModel";
class(cline.meta.rtail.scaled.descending)<-"clineMetaModel";



setupMoleCenterClineParameters<-function(myModel,scaling,x=NULL,y=NULL) {
  pTnames<-names(myModel$parameterTypes);
  if("pMin" %in% pTnames) mdlMin <- "pMin";
  if("xMin" %in% pTnames) mdlMin <- "xMin";
  if("pMax" %in% pTnames) mdlMax <- "pMax";
  if("xMax" %in% pTnames) mdlMax <- "xMax";
  
  if(scaling=="none"){
    attr(myModel$parameterTypes[[mdlMin]],"fixed")<-TRUE;
    attr(myModel$parameterTypes[[mdlMax]],"fixed")<-TRUE;
    myModel$parameterTypes[[mdlMin]]$val<-0;
    myModel$parameterTypes[[mdlMax]]$val<-1;
    
  } else if(scaling=="fixed") {
    attr(myModel$parameterTypes[[mdlMin]],"fixed")<-TRUE;
    attr(myModel$parameterTypes[[mdlMax]],"fixed")<-TRUE;
    if(!is.null(y)){
      myModel$parameterTypes[[mdlMin]]$val<-min(y);
      myModel$parameterTypes[[mdlMax]]$val<-max(y);
    }
  } else if(scaling=="free") {
    attr(myModel$parameterTypes[[mdlMin]],"fixed")<-FALSE;
    attr(myModel$parameterTypes[[mdlMax]],"fixed")<-FALSE;
    if(!is.null(y)){
      myModel$parameterTypes[[mdlMin]]$val<-min(y);
      myModel$parameterTypes[[mdlMax]]$val<-max(y);
      cline.suggestionFunc1D[[mdlMin]](x,y)->junk;
    attr(myModel$parameterTypes[[mdlMin]],"limit.lower")<-junk[[1]];
    attr(myModel$parameterTypes[[mdlMin]],"limit.upper")<-junk[[2]];
      cline.suggestionFunc1D[[mdlMax]](x,y)->junk;
    attr(myModel$parameterTypes[[mdlMax]],"limit.lower")<-junk[[1]];
    attr(myModel$parameterTypes[[mdlMax]],"limit.upper")<-junk[[2]];
    }
  } else {
    stop(paste("Scaling type",scaling,"unrecignized. Please use none, fixed, or free."));
  }
   
  
  if(!is.null(x)){
    qX<-quantile(x,probs=c(0.25,0.5,0.75));
    myModel$parameterTypes$center$val<-qX[[2]]; 
    myModel$parameterTypes$width$val<-qX[[3]]-qX[[1]];
      cline.suggestionFunc1D$center(x,y)->junk;
    attr(myModel$parameterTypes$center,"limit.lower")<-junk[[1]];
    attr(myModel$parameterTypes$center,"limit.upper")<-junk[[2]];
    
      cline.suggestionFunc1D$width(x,y)->junk;
    attr(myModel$parameterTypes$width,"limit.lower")<-junk[[1]];
    attr(myModel$parameterTypes$width,"limit.upper")<-junk[[2]];
    index<-"deltaR";
    if(index %in% pTnames){
      cline.suggestionFunc1D[[index]](x,y)->junk;
      attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
      attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
    }
    index<-"deltaM";
    if(index %in% pTnames){
      cline.suggestionFunc1D[[index]](x,y)->junk;
      attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
      attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
    }
    index<-"deltaL";
    if(index %in% pTnames){
      cline.suggestionFunc1D[[index]](x,y)->junk;
      attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
      attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
    }
    
    
  }
   index<-"tauR";
   if(index %in% pTnames){
     cline.suggestionFunc1D[[index]](x,y)->junk;
     attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
     attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
   }
   index<-"tauM";
   if(index %in% pTnames){
     cline.suggestionFunc1D[[index]](x,y)->junk;
     attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
     attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
   }
   index<-"tauL";
   if(index %in% pTnames){
     cline.suggestionFunc1D[[index]](x,y)->junk;
     attr(myModel$parameterTypes[[index]],"limit.lower")<-junk[[1]];
     attr(myModel$parameterTypes[[index]],"limit.upper")<-junk[[2]];
   }
    
  
  return(myModel);
}

buildCline1D<-function(data,scaling,direction,ascending.cline,descending.cline){
dist=NULL;
obs=NULL;
  if(is.null(direction) ){
    if(is.null(data))
      stop("Either provide observation data, or set direction.");
    if(inherits(data,"clineSampleData1D")){
      if(mean(data$frame$dist) <
         (sum(data$frame$dist*data$frame$obsFreq)/sum(data$frame$obsFreq))){
        direction="ascending";
      } else {
        direction="descending";
      }
      dist=data$frame$dist;
      obs=data$frame$obsFreq;
    } else if(inherits(data,"clineSampleData1DCLT")){
      if(mean(data$frame$dist) <
         (sum(data$frame$dist*data$frame$obsMean)/sum(data$frame$obsMean))){
        direction="ascending";
      } else {
        direction="descending";
      }
      dist=data$frame$dist;
      obs=data$frame$obsMean;
    } else {
      stop("Class of data unknown.");
    }
  }else{
    if(inherits(data,"clineSampleData1D")){
      
      dist=data$frame$dist;
      obs=data$frame$obsFreq;
    } else if(inherits(data,"clineSampleData1DCLT")){
      
      dist=data$frame$dist;
      obs=data$frame$obsMean;
    }
  }
  if(identical(tolower(direction),"ascending")){
    myModel=ascending.cline;
  }else if(identical(tolower(direction),"descending")){
    myModel=descending.cline;
  }else {
    stop(paste("Direction",direction,"not recognized. Specify either ascending or descending."));
  }
  myModel=setupMoleCenterClineParameters(myModel,scaling,dist,obs);
  attr(myModel,"tails")<-"none";
  return(myModel);
}

## makeSimpleCline1D<-function(data=NULL,scaling="none",direction=NULL){
##   return(buildCline1D(data,scaling,direction,
##                       cline.meta.simple.scaled.ascending,
##                       cline.meta.simple.scaled.descending));
## }

## makeTailedCline1D<-function(data=NULL,scaling="none",direction=NULL){
##   myTailedModel<-buildCline1D(data,scaling,direction,
##                               cline.meta.tailed.scaled.ascending,
##                               cline.meta.tailed.scaled.descending);
##   Width<-myTailedModel$parameterTypes$width$val;
##   attr(myTailedModel,"tails")<-"both";
## ##   myTailedModel$parameterTypes$deltaL$val<-Width;
## ##   myTailedModel$parameterTypes$deltaR$val<-Width;
##   return(myTailedModel);
## }

hzar.makeCline1DFreq<- function(data=NULL,scaling="none",tails="none",direction=NULL){
  if(identical(tolower(tails),"none")){
  ##   return(makeSimpleCline1D(data,scaling,direction));
 return(buildCline1D(data,scaling,direction,
                      cline.meta.simple.scaled.ascending,
                      cline.meta.simple.scaled.descending));
  }else if(identical(tolower(tails),"both")) {
##     return(makeTailedCline1D(data,scaling,direction));
      myTailedModel<-buildCline1D(data,scaling,direction,
                              cline.meta.tailed.scaled.ascending,
                              cline.meta.tailed.scaled.descending);
  Width<-myTailedModel$parameterTypes$width$val;
  attr(myTailedModel,"tails")<-"both";
##   myTailedModel$parameterTypes$deltaL$val<-Width;
##   myTailedModel$parameterTypes$deltaR$val<-Width;
  return(myTailedModel);
  }else if(identical(tolower(tails),"right")) {
myRightCline<-buildCline1D(data,scaling,direction,
                              cline.meta.rtail.scaled.ascending,
                              cline.meta.rtail.scaled.descending);
  attr(myRightCline,"tails")<-"right";
    return(myRightCline);
  }else if(identical(tolower(tails),"left")) {
myLeftCline<-buildCline1D(data,scaling,direction,
                              cline.meta.ltail.scaled.ascending,
                              cline.meta.ltail.scaled.descending);
  attr(myLeftCline,"tails")<-"left";
    return(myLeftCline);
  }else if(identical(tolower(tails),"mirror")) {
myMirrorCline<-buildCline1D(data,scaling,direction,
                              cline.meta.mtail.scaled.ascending,
                              cline.meta.mtail.scaled.descending);
  attr(myMirrorCline,"tails")<-"mirror";
    return(myMirrorCline);
  }
  stop(paste("Cline with",tails,"tail(s) not available.")); 
}




## ## format of the cline Multi Model Frame
## objCMMFrame<-list();
## objCMMFrame$modelFrames <- list(objCFrame1,objCFrame2);
## objCMMFrame$data <-objSampleData;  ##Keep Seperate?
## objCMMFrame$parameters <-objParameterDataFrame;  ##Attach to frames?
## class(objCMMFrame)<-"clineMultiModelFrame";
