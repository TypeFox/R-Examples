
model.addReqClause <- function(meta.model,newClause){
  oldClause <- body(meta.model$req)[[2]][[2]];
  body(meta.model$req)[[2]][[2]]<-substitute(clause1 & clause2,list(clause1=oldClause,clause2=newClause));
  return(meta.model);
}

hzar.model.addMaxWidth <- function(meta.model,maxValue){
  clause <- substitute(width<V,list(V=maxValue));
  attr(meta.model$parameterTypes$width,"limit.upper") <- maxValue;
  return(model.addReqClause(meta.model,clause));
}

hzar.model.addMaxCenter <- function(meta.model,maxValue){
  clause <- substitute(center<V,list(V=maxValue));
  attr(meta.model$parameterTypes$center,"limit.upper") <- maxValue;
  return(model.addReqClause(meta.model,clause));
}

hzar.model.addMaxDelta <- function(meta.model,maxValue){
  if(identical(attr(meta.model,"tails"),"none"))
    return(meta.model);
  
  clauseL <- substitute(deltaL<V,list(V=maxValue));
  clauseM <- substitute(deltaM<V,list(V=maxValue));
  clauseR <- substitute(deltaR<V,list(V=maxValue));
  assignLimit <- function(param){
    attr(meta.model$parameterTypes[[param]],"limit.upper") <<- maxValue;
  }
  if(identical(attr(meta.model,"tails"),"both")){
    assignLimit("deltaL");
    assignLimit("deltaR");
    return(model.addReqClause(model.addReqClause(meta.model,clauseL),clauseR));
  }
  if(identical(attr(meta.model,"tails"),"right")){
    assignLimit("deltaR");
    return(model.addReqClause(meta.model,clauseR));
  }
  if(identical(attr(meta.model,"tails"),"left")){
    assignLimit("deltaL");
    return(model.addReqClause(meta.model,clauseL));
  }
  
  if(identical(attr(meta.model,"tails"),"mirror")){
    assignLimit("deltaM");
    return(model.addReqClause(meta.model,clauseM));
  }
  
}
hzar.model.addMinCenter <- function(meta.model,minValue){
  clause <- substitute(center>V,list(V=minValue));
  attr(meta.model$parameterTypes$center,"limit.lower") <- minValue;
  return(model.addReqClause(meta.model,clause));
}

hzar.model.addCenterRange <- function(meta.model,low,high){
  return(hzar.model.addMaxCenter(hzar.model.addMinCenter(meta.model,low),high))
}

hzar.model.addBoxReq <- function(meta.model,low,high){
  return(hzar.model.addCenterRange(hzar.model.addMaxDelta(hzar.model.addMaxWidth(meta.model,high-low),high-low),low,high));
}
model.addMinParam <- function(meta.model,paramQ,paramN,minValue){
  clause <- substitute(P>V,list(P=paramQ,V=minValue));
  attr(meta.model$parameterTypes[[paramN]],"limit.lower") <- minValue;
  return(model.addReqClause(meta.model,clause));
}
model.addMaxParam <- function(meta.model,paramQ,paramN,maxValue){
  clause <- substitute(P<V,list(P=paramQ,V=maxValue));
  attr(meta.model$parameterTypes[[paramN]],"limit.upper") <- maxValue;
  return(model.addReqClause(meta.model,clause));
}

model.addParamRange <-  function(meta.model,paramQ,paramN,minValue,maxValue){
  clauseA <- substitute(P>V,list(P=paramQ,V=minValue));
  attr(meta.model$parameterTypes[[paramN]],"limit.lower") <- minValue;
  clauseB <- substitute(P<V,list(P=paramQ,V=maxValue));
  attr(meta.model$parameterTypes[[paramN]],"limit.upper") <- maxValue;
  return(model.addReqClause(model.addReqClause(meta.model,clauseA),clauseB));
}

hzar.model.addMaxVariance <- function(meta.model,maxValue){
  meta.model <- model.addMaxParam(meta.model,quote(varL),"varL",maxValue)
  meta.model <- model.addMaxParam(meta.model,quote(varH),"varH",maxValue)
  meta.model <- model.addMaxParam(meta.model,quote(varR),"varR",maxValue)
  return(meta.model);
}
hzar.model.addMuRange <- function(meta.model,low,high){
  meta.model <- model.addParamRange(meta.model,quote(muL),"muL",low,high)
  meta.model <- model.addParamRange(meta.model,quote(muR),"muR",low,high)
  return(meta.model);
}
hzar.model.addNormalBox <- function(meta.model,left,right,bottom,top){
  meta.model <- hzar.model.addBoxReq(meta.model,left,right)
  meta.model <- hzar.model.addMuRange(meta.model,bottom,top)
  meta.model <- hzar.model.addMaxVariance(meta.model,(top-bottom)^2)
  return(meta.model);
}
  
