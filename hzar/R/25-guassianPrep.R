obs.g.summary <- function(
                          distance,
                          muObs,
                          varObs,
                          nEff,
                          siteID=paste("P",1:length(distance),sep=""),
                          ylim=NULL){
  list()->obj;
  obj$frame <- na.omit(data.frame(
         site=siteID,
         dist=as.numeric(distance),
         mu=as.numeric(muObs),
         var=as.numeric(varObs),
         nEff=as.numeric(nEff),
         row.names=siteID))
  obj$model.LL <- compileLLfunc(obj)
  if(identical(is.null(ylim),TRUE)){
    #attach(obj$frame)
    f=obj$frame;
    ylim=extendrange(c(min(f$mu-sqrt(f$var/f$nEff)),
      max(f$mu+sqrt(f$var/f$nEff)) ));
    #detach();
  }
  obj$ylim <- ylim;
  class(obj)<-c("guassSampleData1D","hzar.obsData");
  obj
}
hzar.doNormalData1DPops <- obs.g.summary

hzar.mapSiteDist <- function(siteID ,distance){
  distance <- as.numeric(distance);
  names(distance) <- siteID;
  return(distance);
}

hzar.doNormalData1DRaw <- function(site.dist,traitSite,traitValue){
  index <- names(site.dist );
  tV <- lapply(index,function(x) na.omit(traitValue[traitSite==x]))
  ## tmp <- (sapply(tV,length)>0)
  ## index <- index[tmp]; tV <- tV[tmp];site.dist <- site.dist[tmp]
  obs.g.summary(siteID=index,
              distance=site.dist,
              muObs=sapply(tV,mean),
              varObs=sapply(tV,var),
              nEff=sapply(tV,length),
                ylim=extendrange(c(min(na.omit(traitValue)),
                       max(na.omit(traitValue))))
              )
}

compileLLfunc <- function(obsData,
                          muExp=quote(clineMean(x)),
                          varExp=quote(clineVar(x)),
                          baseFunc=function(clineMean,clineVar)0){
 body(baseFunc) <-
   guassianThetaLLExpF(distance=obsData$frame$dist,
                       sampleMean=obsData$frame$mu,
                       sampleVar=obsData$frame$var,
                       nEff=obsData$frame$nEff,
                       muExp=muExp,
                       varExp=varExp)
 environment(baseFunc) <- .GlobalEnv
 baseFunc
}
gCA <- list();
gCA$mu <- alist(muL=,muR=)
gCA$var <- alist(varL=,varR=,varH=)
gCA$none <- c(gCA$mu,gCA$var,alist(center=,width=))
gCA$mirror <- c(gCA$none,alist(deltaM=,tauM=))
gCA$right <- c(gCA$none,alist(deltaR=,tauR=))
gCA$left <- c(gCA$none,alist(deltaL=,tauL=))
gCA$both <- c(gCA$none,alist(deltaL=,tauL=,deltaR=,tauR=))

makeGExp <- function(muExp,vExp) list(muExp=muExp,vExp=vExp)
step1VGExpF <-  function( cExp, trueGExp, falseGExp )
  makeGExp(step1VectorExpF(cExp,trueGExp$muExp,falseGExp$muExp),
           step1VectorExpF(cExp,trueGExp$vExp,falseGExp$vExp))
mvGExpF <- function( clineExp, vClineExp,
                    swap=FALSE,
                    minMu=quote(muL),
                    maxMu=quote(muR),
                    leftVar=quote(varL),
                    rightVar=quote(varR),
                    kappaE=quote(varH)){
  if(swap){
    return(makeGExp(muExpF(maxMu,minMu,clineExp),
                    varExpF(rightVar,leftVar,kappaE,clineExp,vClineExp)))
  }else{
    return(makeGExp(muExpF(minMu,maxMu,clineExp),
                    varExpF(leftVar,rightVar,kappaE,clineExp,vClineExp)
                    ))
  }
}

gTailExpF <- function(LdeltaE,tauE,reverse=FALSE,swap=reverse){
  if(reverse){
    return(mvGExpF(pTailExpF(LdeltaE,tauE,clineLogitRev),
                   vTailExpF(LdeltaE,tauE,clineLogitRev),
                   swap=swap))
  }else{
    return(mvGExpF(pTailExpF(LdeltaE,tauE,clineLogit),
                   vTailExpF(LdeltaE,tauE,clineLogit),
                   swap=swap))
  }
}
gCenterExpF <- function(reverse=FALSE,swap=FALSE){
  if(reverse){
    return(mvGExpF(sigmoidExpF(clineLogit),
                   vCenterExpF(clineLogit),
                   swap=swap))
  }else{
    return(mvGExpF(sigmoidExpF(clineLogitRev),
                   vCenterExpF(clineLogit),
                   swap=swap))
  }
}


g.suggest <- list();
g.suggest$muL <- function(data,...) data$mu[which.min(data$dist)]
g.suggest$muR <- function(data,...) data$mu[which.max(data$dist)]
g.suggest$varL <- function(data,...) {
  data2 <- data[data$var > 0.10,c("dist","var")]
  if(nrow(data2)>0) return(max(1,data2$var[which.min(data2$dist)]));
  return(1);
}
g.suggest$varR <- function(data,...){
  data2 <- data[data$var > 0.10,c("dist","var")]
  if(nrow(data2)>0) return(max(1,data2$var[which.max(data2$dist)]));
  return(1);
}
g.suggest$varH <- function(data,varL,varR,...) max(1,max(data$var)-(varR+varL)/2)
g.suggest$center <- function(data,...) data$dist[which.max(data$var)]
g.suggest$width <- function(data,muL,muR,...) {
  dist <- data$dist
  res <- -1
  if(muL<muR){
    res <- (min(dist[(5* data$mu >= (4*muR+muL))])
           -max(dist[(5* data$mu <= (4*muL+muR))]))
    
  }else if(muL>muR){
    res <- (min(dist[(5* data$mu <= (4*muR+muL))])
           -max(dist[(5* data$mu >= (4*muL+muR))]))
                 

  }
  if(all(!is.na(res),!is.nan(res),res>0))return(res)
  if(length(dist)==0)return(1)
  dist <- sort(unique(dist))
  l <- length(dist)
  if(l==1)return(dist)
  if(l<4)return(max(dist)-min(dist)+1)  
  return(min(dist[3:l]-dist[1:(l-2)]))
}
g.suggest$deltaL <- function(width,...)3*width/4
g.suggest$deltaM <- function(width,...)3*width/4
g.suggest$deltaR <- function(width,...)3*width/4
g.suggest$tauL <- function(...) 0.5
g.suggest$tauM <- function(...) 0.5
g.suggest$tauR <- function(...) 0.5

g.suggestAll <- function(obsData,mArgs){
  index <- intersect(names(g.suggest),names(mArgs))
  dataL <- alist(data=obsData$frame)
  for(iter in index)
    mArgs[[iter]] = do.call(g.suggest[[iter]],c(dataL,mArgs))
  mArgs
}

g.sMuVarL <- function(obsData){
  data <- obsData$frame
  data2 <- data[data$var > 0,c("dist","var","nEff")]
  if(nrow(data2)>0) {
    i <- which.min(data2$dist)
    return(data2$var[i]/data2$nEff[i]);
  } else {
    return(var(data$mu)/sum(data$nEff));
  }
}

g.sMuVarR <- function(obsData){
  data <- obsData$frame
  data2 <- data[data$var > 0,c("dist","var","nEff")]
  if(nrow(data2)>0) {
    i <- which.max(data2$dist)
    return(data2$var[i]/data2$nEff[i]);
  } else {
    return(var(data$mu)/sum(data$nEff));
  }
}

hzar.makeCline1DNormal <- function(data,tails="none"){
  model <- list(mu=0,
                var=1,
                
                args=alist(),
                req=function(varL,varR,varH,width)
                  return((width>0)&(varL>0)&(varR>0)&(varH>0)),
                mFunc=function()0,
                vFunc=function()0,
                init=list())
  class(model) <- "clineMetaModel";
  
  if(identical(tolower(tails),"none")){
    mV <- gCenterExpF();
    model$mu=mV$muExp
    model$var=mV$vExp
    model$args=gCA$none
    attr(model,"tails")<-"none";
  }else if(identical(tolower(tails),"left") ){
    mV <- step1VGExpF(quote(x < center  - deltaL),
                      gTailExpF(quote(4*deltaL/width),quote(tauL)),
                      gCenterExpF());
    model$mu=mV$muExp
    model$var=mV$vExp
    model$args=gCA$left
    model=model.addReqClause(model,
      quote((deltaL>0)&(tauL>=0)&(tauL<=1)))
    attr(model,"tails")<-"left";
  }else if(identical(tolower(tails),"right")){
    mV <- step1VGExpF(quote(x > center  + deltaR),
                      gTailExpF(quote(4*deltaR/width),quote(tauR),TRUE),
                      gCenterExpF());
    model$mu=mV$muExp
    model$var=mV$vExp
    model$args=gCA$right
    model=model.addReqClause(model,
      quote((deltaR>0)&(tauR>=0)&(tauR<=1)))
    attr(model,"tails")<-"right";

  }else if(identical(tolower(tails),"both") ){
mV <- step1VGExpF(quote(x < center  - deltaL),
                  gTailExpF(quote(4*deltaL/width),quote(tauL)),
                  step1VGExpF(quote(x > center  + deltaR),
                              gTailExpF(quote(4*deltaR/width),
                                        quote(tauR),
                                        TRUE),
                              gCenterExpF()));
    model$mu=mV$muExp
    model$var=mV$vExp
    model$args=gCA$both
    model=model.addReqClause(model,
      quote((deltaL>0)&(tauL>=0)&(tauL<=1)))
    model=model.addReqClause(model,
      quote((deltaR>0)&(tauR>=0)&(tauR<=1)))
    attr(model,"tails")<-"both";

  }else if(identical(tolower(tails),"mirror")){
mV <- step1VGExpF(quote(x < center  - deltaM),
                  gTailExpF(quote(4*deltaM/width),quote(tauM)),
                  step1VGExpF(quote(x > center  + deltaM),
                              gTailExpF(quote(4*deltaM/width),
                                        quote(tauM),
                                        TRUE),
                              gCenterExpF()));
    model$mu=mV$muExp
    model$var=mV$vExp
    model$args=gCA$mirror
    
    model=model.addReqClause(model,
      quote((deltaM>0)&(tauM>=0)&(tauM<=1)))
    attr(model,"tails")<-"mirror";

  }else {
    stop("Unrecognized type requested.");
  }
  model$parameterTypes <- CLINEPARAMETERS[intersect(names(CLINEPARAMETERS),
                                                    names(model$args))]
  meta.init(model) <- g.suggestAll(data,model$args)
  meta.fix(model) <- FALSE;#rep(FALSE,length(model$init))
  ## names(model$fixed) <- names( model$init)
  meta.tune(model) <- 1.5; #$tune <- as.list(rep(1.5,length(model$init)))
  ## names(model$tune) <- names( model$init)
  
  formals(model$req,envir=.GlobalEnv) <- model$args
  formals(model$mFunc,envir=.GlobalEnv) <- model$args
  formals(model$vFunc,envir=.GlobalEnv) <- model$args
  mS <- substitute(body(res) <- substitute(mu),model)
  vS <- substitute(body(res) <- substitute(var),model)
  body(model$mFunc) <- as.call(c(as.name("{"),
                                 expression(res <- function(x) x),
                                 mS,expression(environment(res) <- .GlobalEnv,
                                     return(res))))
  body(model$vFunc) <- as.call(c(as.name("{"),
                                 expression(res <- function(x) x),
                                 vS,expression(environment(res) <- .GlobalEnv,
                                     return(res))))
  muVarL <- g.sMuVarL(data);
  muVarR <- g.sMuVarR(data);
  meta.lower(model)$muL <- data$ylim[[1]]
  meta.lower(model)$muR <- data$ylim[[1]]
  meta.upper(model)$muL <- data$ylim[[2]]
  meta.upper(model)$muR <- data$ylim[[2]]
##   attr(model$parameterTypes$muL,"limit.lower") <- data$ylim[[1]]#model$init$muL-sqrt(muVarL)
##   attr(model$parameterTypes$muL,"limit.upper") <- data$ylim[[2]]#model$init$muL+sqrt(muVarL)
##   attr(model$parameterTypes$muR,"limit.lower") <- data$ylim[[1]]#model$init$muR-sqrt(muVarR)
##   attr(model$parameterTypes$muR,"limit.upper") <- data$ylim[[2]]#model$init$muR+sqrt(muVarR)
  meta.lower(model)$varL<- min(model$init$varL/100,1e-8)
  meta.upper(model)$varL<- max(model$init$varL*1e4,1e8)
  meta.lower(model)$varH<- min(model$init$varH/100,1e-8)# model$init$varH/100
  meta.upper(model)$varH<- max(model$init$varH*1e4,1e8)
  meta.lower(model)$varR<-  min(model$init$varR/100,1e-8)#model$init$varR/100
  meta.upper(model)$varR<- max(model$init$varR*1e4,1e8)
  model
}
