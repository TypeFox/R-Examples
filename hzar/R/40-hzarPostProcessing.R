## So, useful tests and helper functions:

## Are these two fit requests from the same cline model?
hzar.sameModel <- function(fitA, fitB){
  if(!inherits(fitA, c("hzar.fitRequest","hzar.dataGroup")) ||
     !inherits(fitB, c("hzar.fitRequest","hzar.dataGroup")) ){
    stop("Can only compare hzar.fitRequest or hzar.dataGroup objects.");
  }
  if(!hzar.sameObsData(fitA,fitB))
    return(FALSE);
  chkA=cline.extract.modelPrep(fitA);
  chkB=cline.extract.modelPrep(fitB);
  return(all(identical(body(chkA$model.gen),body(chkB$model.gen)),
             identical(chkA$new.formals,chkB$new.formals) ));

}

## Are these two objects from the same observation data?
hzar.sameObsData <- function(fitA, fitB){
  if(!inherits(fitA, c("hzar.fitRequest",
                       "hzar.dataGroup",
                       "clineSampleData1D",
                       "clineSampleData1DCLT",
                       "hzar.obsData",
                       "hzar.obsDataGroup")) ||
     !inherits(fitB, c("hzar.fitRequest",
                       "hzar.dataGroup",
                       "clineSampleData1D",
                       "clineSampleData1DCLT",
                       "hzar.obsData",
                       "hzar.obsDataGroup")) ){
    stop("Can only compare hzar data objects.");
  }
  return(identical(hzar.extract.obsData(fitA)$frame,
                   hzar.extract.obsData(fitB)$frame ));

}

hzar.copyModelLabels<-function(group1,group2){
  if(inherits(group1,c("hzar.obsDataGroup")))
    group1<-group1$data.groups;
  if(inherits(group2,c("hzar.obsDataGroup"))){
    group2$data.groups<-hzar.copyModelLabels(group1,group2$data.groups);
    return(group2);
  }
  if(!is.list(group1)||!is.list(group2))
    stop("group1 or group2 of malformed type");
  
  names(group2)<-
    as.character(lapply(group2,
                        function(x)
                        names(group1)
                        [which(as.logical(lapply(group1,
                                                 hzar.sameModel,
                                                 x)))[1]]));
  return(group2);
}
  

## hzar.modelLabel.set

## works on: clineMetaModel, hzar.fitRequest, hzar.dataGroup,
## hzar.cline, hzar.obsDataGroup[multiple unique].

## hzar.modelLabel.get                     

## hzar.traitLabel.set                     

## works on: hzar.obsData, hzar.fitRequest, hzar.dataGroup,
## hzar.obsDataGroup, hzar.cline.

## hzar.traitLabel.get                     


## Get free parameter count (for several model object types?)

## Extract meta model object?
## as.list(environment(fitRequest$llFunc))

cline.extract.modelPrep <- function(fitRequest){
  if(inherits(fitRequest, c("hzar.fitRequest","hzar.dataGroup"))){
    llFunc<-fitRequest$llFunc;
  }else if(is.function(fitRequest)){
    llFunc<-fitRequest;
  } else{
    stop(paste("hzar does not understand objects of type ",class(fitRequest)));
  }
  return(as.list(environment(llFunc)));
}
hzar.extract.old.model.gen <- function(fitRequest){
  ##returns a function of TTT, the list of free parameters
  
  return(cline.extract.modelPrep(fitRequest)$model.gen);
}

hzar.extract.old.model.req <- function(fitRequest){
  ##returns a function of TTT, the list of free parameters
  return(cline.extract.modelPrep(fitRequest)$model.req);
}

## hzar.extract.old.model.prior <- function(fitRequest){
##   ##returns a function of TTT, the list of free parameters
## }

## Extract cline distance function?

## Extract observed data

hzar.extract.obsData <- function(fitRequest) {
  if(inherits(fitRequest,c("hzar.dataGroup","hzar.obsDataGroup"))){
    return(fitRequest$obsData);
  }
  if(inherits(fitRequest,c("clineSampleData1D","clineSampleData1DCLT","hzar.obsData"))){
    return(fitRequest);
  }
      
  return(cline.extract.modelPrep(fitRequest)$obsData);
}

## Cline Object (free parameter values, parameter values, cline
## function of distance, likelihood value)

hzar.make.cline<-function(free.parameters,parameters,func,LL,isValid=is.function(func)){
  ##simple packaging
  obj<-list(param.free=free.parameters,param.all=parameters,
            clineFunc=func, logLike=LL,isValid=isValid);
  class(obj)<-"hzar.cline";
  return(obj);
}

## generate hzar.cline object for parameters in the given context

hzar.gen.cline<-function(free.parameters,fitRequest){
  ## print(free.parameters);
  context<-cline.extract.modelPrep(fitRequest);
  cline.func<-NULL;
  cline.param<-c(as.list(free.parameters),context$param.fixed)
  if(!all(identical(formals(context$model.req),context$new.formals),
          identical(formals(context$model.gen),context$new.formals))){
    formals(context$model.req) <- context$new.formals
    formals(context$model.gen) <- context$new.formals
  }
     
  if(do.call(context$model.req,as.list(free.parameters))){
    cline.func<-do.call(context$model.gen,as.list(free.parameters));
  }
  return(hzar.make.cline(free.parameters,parameters=cline.param,func=cline.func,LL=context$llFunc(free.parameters)));
}
## Extract ML cline object

hzar.get.ML.cline <- function(fitRequest){
  ## check for object type?
  if( inherits(fitRequest,"hzar.fitRequest") ){
    ## check for evaluation (aka fitting)
    if(! identical( attr(fitRequest,"fit.success") , TRUE)){
      stop("fitRequest needs a successful MCMC sampling chain.");
    }
    ## get list of logLikelihoods
    data.LL<-hzar.eval.clineLL(llFunc=fitRequest$llFunc,
                               data=fitRequest$mcmcRaw);
    data.param<-as.data.frame(fitRequest$mcmcRaw);
  } else if(inherits(fitRequest,"hzar.dataGroup")){
    if(!identical(is.null(fitRequest$ML.cline),TRUE))
      return(fitRequest$ML.cline);
    data.LL<-fitRequest$data.LL;
    data.param<-fitRequest$data.param;
  } else {
    stop(paste("hzar does not understand objects of type ",
               class(fitRequest)));
  }
  ## get parameters corresponding to max value
  param.free = data.param[data.LL==max(data.LL), ,drop=FALSE][1, ,drop=FALSE];
  ## generate cline
  return(hzar.gen.cline(param.free,fitRequest));
}

## convenience classes to gather post-processing data

## process mcmc data from the same model and the same obsData

hzar.make.dataGroup<-
  function(data.mcmc,llFunc, ML.cline=NULL, doPar=FALSE,
           data.LL=hzar.eval.clineLL(llFunc=llFunc,data=data.mcmc,doPar=doPar),
           data.param=as.data.frame(data.mcmc),
           obsData=hzar.extract.obsData(llFunc)){
  result<-list(llFunc=llFunc,
               data.mcmc=data.mcmc,
               data.param=data.param,
               data.LL=data.LL,
               ML.cline=ML.cline,
               obsData=obsData);
  class(result)<-"hzar.dataGroup";
  result$ML.cline<-hzar.get.ML.cline(result);
  return(result);
}

hzar.fit2DataGroup<-function(fitRequest,doPar=FALSE){
  ## check for object type?
  if( inherits(fitRequest,"hzar.fitRequest") ){
    ## check for evaluation (aka fitting)
    if(! identical( attr(fitRequest,"fit.success") , TRUE)){
      ## stop("fitRequest needs a successful MCMC sampling chain.");
      mcmcPsuedo<- mcmc(as.data.frame(fitRequest$modelParam$init));
      return(hzar.make.dataGroup(data.mcmc=mcmcPsuedo,
                                 llFunc<-fitRequest$llFunc));
    }
    return(hzar.make.dataGroup(data.mcmc=fitRequest$mcmcRaw,
                               llFunc<-fitRequest$llFunc,
                               doPar=doPar));
  } else if(inherits(fitRequest,"hzar.dataGroup")){
    return(fitRequest);
  } else {
    stop(paste("hzar does not understand objects of type ",
               class(fitRequest)));
  }
}

## collating function

## need to add version that works on a list of runs with no other
## arguments.  I should also consider adding a class the describes a
## set of chains of sequential runs to the fitting methods.

hzar.dataGroup.add<-function(dataGroup, fitRequestL=list(),doPar=FALSE){
  if( inherits(dataGroup, c("hzar.fitRequest","hzar.dataGroup")) ){
    dataGroup<-hzar.fit2DataGroup(dataGroup,doPar=doPar);
  } else if(is.list(dataGroup) && (length(fitRequestL)==0)){
    if(length(dataGroup)==0)
      return(NULL);
    oldFitRequest<-dataGroup;
    dataGroup<-hzar.fit2DataGroup(dataGroup[[1]],doPar=doPar);
    if(length(oldFitRequest)>1)
      fitRequestL<-lapply(2:length(oldFitRequest),
                         function(x)oldFitRequest[[x]]);
  }else {
    stop("dataGroup not of apropriate type.");
  }
  if( inherits(fitRequestL, c("hzar.fitRequest","hzar.dataGroup")) ){
    fitRequestL<-hzar.fit2DataGroup(fitRequestL,doPar=doPar);
  } else if( is.list(fitRequestL) ){
    lapply(fitRequestL,function(x)
           dataGroup<<-hzar.dataGroup.add(dataGroup,x,doPar=doPar));
    return(dataGroup);
  } else {
    stop("fitRequestL not of appropriate type.");
  }
  if( !hzar.sameModel(dataGroup,fitRequestL))
    stop("dataGroup and fitRequestL must use the same model.");
  Left.Cline<-hzar.get.ML.cline(dataGroup);
  Right.Cline<-hzar.get.ML.cline(fitRequestL);
  if(Right.Cline$logLike>Left.Cline$logLike)
    Left.Cline<-Right.Cline;
  return(hzar.make.dataGroup(data.mcmc=rbind(dataGroup$data.mcmc,
                               fitRequestL$data.mcmc),
                             llFunc=dataGroup$llFunc,
                             ML.cline=Left.Cline,
                             data.LL=rbind(dataGroup$data.LL,
                               fitRequestL$data.LL),
                             data.param=rbind(dataGroup$data.param,
                               fitRequestL$data.param),
                             obsData=dataGroup$obsData));
}
## a collect of hzar.dataGroup objects that share the same obsData

hzar.make.obsDataGroup<-function(dataGroups,obsDataGroup=NULL){
  if(identical(is.null(obsDataGroup),TRUE)){
    
    if(inherits(dataGroups,c("hzar.obsDataGroup"))){
      ##print("oDG");
      return(dataGroups);
    }
    ##print("dG");
    if(inherits(dataGroups, c("hzar.fitRequest","hzar.dataGroup"))){
      dataGroups<-list(hzar.fit2DataGroup(dataGroups));
       
      obsData<-hzar.extract.obsData(dataGroups[[1]])
      obsDataGroup<-list(data.groups=dataGroups,
                       obsData=obsData);
      class(obsDataGroup)<-"hzar.obsDataGroup";
      return(obsDataGroup);
    }
    if(inherits(dataGroups,c("hzar.obsData"))){
      obsDataGroup<-list(data.groups=list(),obsData=dataGroups);
      class(obsDataGroup)<-"hzar.obsDataGroup";
      return(obsDataGroup);
    }
    if(is.list(dataGroups)){
      print(length(dataGroups));
      if(inherits(dataGroups[[1]],
                  c("hzar.fitRequest","hzar.dataGroup"))){
        obsDataGroup<-hzar.make.obsDataGroup(dataGroups[[1]]);
      }else if(inherits(dataGroups[[1]],
                        c("hzar.obsDataGroup"))){
        obsDataGroup<-dataGroups[[1]];
      } else {
        return(hzar.make.obsDataGroup(lapply(dataGroups,
                                             hzar.make.obsDataGroup)));
      }
      if(length(dataGroups)==1)
        return(obsDataGroup);
      otherDataGroups<-lapply(2:length(dataGroups),
                              function(x,y) {
                                if(inherits(y[[x]], c("hzar.fitRequest","hzar.dataGroup")))
                                  return(hzar.fit2DataGroup(y[[x]]));
                                y[[x]]
                              },
                              y=dataGroups);
      return(hzar.make.obsDataGroup(otherDataGroups,obsDataGroup));
    }
      
  } else if(inherits(obsDataGroup,c("hzar.obsDataGroup"))){
    if(inherits(dataGroups,c("hzar.obsDataGroup"))){
      return(hzar.make.obsDataGroup(dataGroups$data.groups,obsDataGroup));
    }
    if(inherits(dataGroups, c("hzar.fitRequest","hzar.dataGroup"))){
      dataGroups<-list(hzar.fit2DataGroup(dataGroups));
    } else if(is.list(dataGroups)&&
              inherits(dataGroups[[1]],
                       c("hzar.fitRequest","hzar.dataGroup"))) {
      if(!all(isDG <- sapply(dataGroups,
                             inherits,
                             c("hzar.fitRequest","hzar.dataGroup")))){
        obsDataGroup <-
          hzar.make.obsDataGroup(lapply(which(isDG),
                                        function(x) hzar.fit2DataGroup(dataGroups[[x]])),
                                 obsDataGroup)
        dataGroups <- lapply(which(!isDG),function(x) dataGroups[[x]])
      } else {
        dataGroups <- lapply(dataGroups,hzar.fit2DataGroup)
      }
    }
    if(! hzar.sameObsData(obsDataGroup,dataGroups[[1]]) )
      stop("All dataGroups must be from the same observation data." );
    if(!inherits(dataGroups[[1]], c("hzar.fitRequest","hzar.dataGroup"))){
      if(inherits(dataGroups[[1]],c("hzar.obsDataGroup")))
        return(hzar.make.obsDataGroup(hzar.make.obsDataGroup(dataGroups),
                                      obsDataGroup));
      return(hzar.make.obsDataGroup(lapply(dataGroups,
                                           hzar.make.obsDataGroup),
                                    obsDataGroup));
    }
    modelsKnown<-as.logical(lapply(obsDataGroup$data.groups,
                                   hzar.sameModel,fitB=dataGroups[[1]]));
    if(sum(modelsKnown)==0){
      obsDataGroup$data.groups<-c(obsDataGroup$data.groups,
                                 list(dataGroups[[1]]));
    } else if(sum(modelsKnown)==1){
      index=which(modelsKnown);
      obsDataGroup$data.groups[[index]]<-
        hzar.dataGroup.add(obsDataGroup$data.groups[[index]],
                           dataGroups[[1]]);
    } else {
      stop("hzar.make.obsDataGroup found a fatal error in the obsDataGroup structure.");
    }
    if(length(dataGroups)==1)
        return(obsDataGroup);
    otherDataGroups<-lapply(2:length(dataGroups),
                            function(x,y) y[[x]],
                            y=dataGroups);
    return(hzar.make.obsDataGroup(otherDataGroups,obsDataGroup));
        
  }
  stop("Argument types unhandled by hzar.make.obsDataGroup");
    
}



## obsData<-hzar.extract.obsData(dataGroups[[1]])
##   dataGroups<-lapply(dataGroups, function(x) {
##     if(! hzar.sameObsData(obsData,x) )
##       stop("All dataGroups must be from the same observation data." );
##     hzar.fit2DataGroup(x);
##   });

## model testing / scores

## Likelihood ratio tests?

## AIC

## Simple
hzar.AIC.default <- function(maxLL,param.count){
  return(2*(param.count-maxLL));
}

hzar.AICc.default <- function(maxLL,param.count,nObs){
  return(2*(param.count
            - maxLL
            + (param.count
               * (param.count+1)
               / (nObs-param.count-1))));
}

## hzar.cline

hzar.AIC.hzar.cline <- function(cline){
  return(hzar.AIC.default(cline$logLike,length(cline$param.free)));
}

hzar.AICc.hzar.cline <- function(cline,nObs){
  return(hzar.AICc.default(cline$logLike,length(cline$param.free),nObs));
}

## hzar.dataGroup

hzar.AIC.hzar.dataGroup <- function(dataGroup){
  return(hzar.AIC.hzar.cline(hzar.get.ML.cline(dataGroup)));
}

hzar.AICc.hzar.dataGroup <- function(dataGroup){
  return(hzar.AICc.hzar.cline(hzar.get.ML.cline(dataGroup),
                              sum(dataGroup$obsData$frame$n)));
}

## hzar.obsDataGroup [returns a data frame]

hzar.AIC.hzar.obsDataGroup <- function(obsDataGroup,label="AIC",show.count=FALSE,show.param=FALSE){
  model.clines<-lapply(obsDataGroup$data.groups,hzar.get.ML.cline);
  scores<-as.numeric(lapply(model.clines,hzar.AIC.hzar.cline));

  result<-data.frame(AIC=scores);
  if(!identical(is.null(label),TRUE)){
    names(result)<-label;
    colnames(result)<-label;
  }
  if(!identical(is.null(names(obsDataGroup$data.groups)),TRUE)){
    rownames(result)<-names(obsDataGroup$data.groups);
  }
  if(show.count){
    count=as.numeric(lapply(model.clines,function(x)length(x$param.free)));
    result<-cbind(result,count=count);
  }
  ## if(show.param){
  ##   all.names<-unique(unlist(lapply(model.clines,function(x)names(x$param.all))));
  ##   print(all.names);
  ##   for(junk in 1:length(model.clines)){
  ##     ##print( names(model.clines[[junk]]$param.all));
  ##     print(all.names[all.names %in% names(model.clines[[junk]]$param.all)]);
  ##     print(as.numeric(model.clines[[junk]]$param.all));
  ##     result[junk,all.names[all.names %in% names(model.clines[[junk]]$param.all)] ]<- as.numeric(model.clines[[junk]]$param.all);
  ##     print(as.numeric(model.clines[[junk]]$param.all));
  ##   }
  ## }
    return(result);
}

hzar.AICc.hzar.obsDataGroup <- function(obsDataGroup,label="AICc",show.count=FALSE,show.param=FALSE){
  model.clines<-lapply(obsDataGroup$data.groups,hzar.get.ML.cline);
  scores<-as.numeric(lapply(model.clines,
                            hzar.AICc.hzar.cline,
                            sum(obsDataGroup$obsData$frame$n)));

  result<-data.frame(AIC=scores);
  if(!identical(is.null(label),TRUE)){
    names(result)<-label;
    colnames(result)<-label;
  }
  if(!identical(is.null(names(obsDataGroup$data.groups)),TRUE)){
    rownames(result)<-names(obsDataGroup$data.groups);
  }
  if(show.count){
    count=as.numeric(lapply(model.clines,function(x)length(x$param.free)));
    result<-cbind(result,count=count);
  }
    return(result);
}

## Bayes Factors


## Calculate credible likelihood cut
hzar.getCredCut<-function(dataGroup,rejectionPercent=0.05){
  if(inherits(dataGroup, c("hzar.fitRequest","hzar.dataGroup"))){
    dataGroup<-hzar.fit2DataGroup(dataGroup);
  } else if(is.list(dataGroup)&& length(dataGroup)>0&& 
            inherits(dataGroup[[1]],
                     c("hzar.fitRequest","hzar.dataGroup"))){
    warning("Only calculating LL cut for first element in list");
    dataGroup<-hzar.fit2DataGroup(dataGroup[[1]]);
  }
  model.relLL=exp(sort(dataGroup$data.LL$model.LL-dataGroup$ML.cline$logLike));
  credibleLLspace<-data.frame(LL=sort(dataGroup$data.LL$model.LL),
                              percentile=cumsum(model.relLL/sum(model.relLL)));

  credible.LLcut<-min(subset(credibleLLspace,
                             credibleLLspace$percentile>rejectionPercent)$LL);
  return(credible.LLcut);
}


hzar.getCredParam <- function(dataGroup,rejectionPercent=0.05){
  if(inherits(dataGroup, c("hzar.fitRequest","hzar.dataGroup"))){
    dataGroup<-hzar.fit2DataGroup(dataGroup);
  } else if(is.list(dataGroup)&& length(dataGroup)>0&& 
            inherits(dataGroup[[1]],
                     c("hzar.fitRequest","hzar.dataGroup"))){
    warning("Only calculating LL cut for first element in list");
    dataGroup<-hzar.fit2DataGroup(dataGroup[[1]]);
  }
  credible.LLcut<-hzar.getCredCut(dataGroup,rejectionPercent);
  return(dataGroup$data.param[dataGroup$data.LL$model.LL>credible.LLcut,]);
}


PP.fzCline.reduce <- function(candidates,
                                checkFunc=PP.fz.getCheckFunc(names(candidates)),
                                noReCheck=FALSE,glitz=character(0),wedgeCut=8000*length(names(candidates))){
  candidates<-as.data.frame(candidates);
  dim(candidates)[[1]]->numCandidates;
  if(numCandidates>wedgeCut){
    #cat("/");
    a <- PP.fzCline.reduce(candidates[1:as.integer(numCandidates/2),],
                             checkFunc=checkFunc, noReCheck=noReCheck,
                             glitz=c(glitz," "));
    #cat("V");
    b <- PP.fzCline.reduce(candidates[as.integer(1+(numCandidates/2)):numCandidates,],
                             checkFunc=checkFunc, noReCheck=noReCheck,
                             glitz=c(glitz,"|"));
    #cat(":")#,as.integer(log(
    numA <- dim(a)[[1]]#)/log(2)),"+",as.integer(log(
    numB <- dim(b)[[1]]#)/log(2)),sep="");
    ## if(numA<numB){
    ##     res <- (PP.fzCaMnRc(a,#rbind(a[1:(numA-1),],b[numB,]),
    ##                           b[rev(1:(numB)),],
    ##                           checkFunc=checkFunc));
    ##   }else{
    ##     res <- (PP.fzCaMnRc(b,a[rev(1:(numA)),],
    ##                                     #rbind(a[numA,],b[numB,]),
    ##                                     #rbind(b[rev(1:(numB-1)),],a[rev(1:(numA-1)),]),
    ##                           checkFunc=checkFunc));
    ##   }
      

    
    if((length(glitz)==0)||((rev(glitz)[1])==" ")){
      rN <- 1;
      rD <- 2;
      if(numA<numB){
        cutV <- (rN*numA) %/% (rD);
        res <- (PP.fzCaMnRc(rbind(b[numB,],a[cutV:numA,]),
                              rbind(a[rev(1:(cutV-1)),],b[rev(1:(numB-1)),]),
                              checkFunc=checkFunc));
      }else{
        cutV <- (rN*numB) %/% (rD);
        res <- (PP.fzCaMnRc(rbind(a[numA,],b[cutV:numB,]),
                              rbind(b[rev(1:(cutV-1)),],a[rev(1:(numA-1)),]),
                              checkFunc=checkFunc));
      }
      
      cat("+");
      if(numCandidates>16*wedgeCut)
        cat("\n");
    } else {
      if(numA<numB){
        res <- (PP.fzCaMnRc(a,#rbind(a[1:(numA-1),],b[numB,]),
                              b[rev(1:(numB)),],
                              checkFunc=checkFunc));
      }else{
        res <- (PP.fzCaMnRc(b,a[rev(1:(numA)),],
                                        #rbind(a[numA,],b[numB,]),
                                        #rbind(b[rev(1:(numB-1)),],a[rev(1:(numA-1)),]),
                              checkFunc=checkFunc));
      }
      

      cat("-");
      
    }
    
    return(res);

    
    
  }
  ## cat(glitz,"-",sep="");
##print(numCandidates);
  if(numCandidates<8)
    return(candidates);
  
  return(PP.fzCaMnRc(candidates[1:4,],
                          candidates[5:numCandidates,],
                          checkFunc=checkFunc));
}
 
PP.fzCline.add <-function(accepted,
                            candidate,
                            checkFunc=PP.fz.getCheckFunc(names(accepted)) ){
   if(any(apply(candidate[,names(accepted)],2,max)>apply(accepted,2,max),
          apply(candidate[,names(accepted)],2,min)<apply(accepted,2,min)))
     return(PP.fzCline.addBest (accepted,candidate,checkFunc));
  
  
  res<-PP.fzCline.addMany(accepted=accepted, candidates=candidate,
                            checkFunc=checkFunc);
  return(res);
}
PP.fzCline.addOne <-function(accepted,
                            candidate,
                            checkFirst=1,
                            checkFunc=PP.fz.getCheckFunc(names(accepted)) ){
  accepted<-as.data.frame(accepted);
  candidate<-as.data.frame(candidate);
  dim(accepted)[[1]]->numAccepted;
  ## Check for multiple candidates
  if(dim(candidate)[[1]]>1)
    return(PP.fzCline.addOne(PP.fzCline.addOne(accepted,candidate[1,],
                                                   checkFunc=checkFunc,
                                                   checkFirst=checkFirst),
                               candidate[2:(dim(candidate)[[1]]),],
                               checkFunc=checkFunc,
                               checkFirst=(checkFirst+1)%%4));
  ##rejection?
  if(checkFunc(accepted,candidate))
    accepted<-rbind(accepted,candidate);
  dim(accepted)[[1]]->numAccepted;

  return(accepted);
}

PP.fzCline.addMany <- function(accepted,candidates,checkFunc){
  ## First, check alignment, in order to simplify tail call
  numCan <- dim(candidates)[[1]]
  if(numCan<4)
    return(PP.fzCline.addOne(accepted,
                               candidate=candidates,
                               checkFunc=checkFunc));
 
  ## Assumes all intelligent work has been done, so brute force.
  ## Unroll loops, too.
  for(canI in 4*(1:(numCan%/%4))-3){
    if(checkFunc(accepted,candidates[canI,]))
      accepted<-rbind(accepted,candidates[canI,]);
    if(checkFunc(accepted,candidates[canI+1,]))
      accepted<-rbind(accepted,candidates[canI+1,]);
    if(checkFunc(accepted,candidates[canI+2,]))
      accepted<-rbind(accepted,candidates[canI+2,]);
    if(checkFunc(accepted,candidates[canI+3,]))
      accepted<-rbind(accepted,candidates[canI+3,]);
    if((dim(accepted)[[1]]->numAccepted)>16){
      if(checkFunc(accepted[2:numAccepted,],accepted[1,])){
        accepted <- accepted[c(2:numAccepted,1),];
      }else{
        accepted <- accepted[2:numAccepted,];
      }
    }    
  }
  canI= 4*(numCan%/%4)-3
  if((!numCan<(canI+4))&&checkFunc(accepted,candidates[canI+4,]))
    accepted<-rbind(accepted,candidates[canI+4,]);
  if((!numCan<(canI+5))&&checkFunc(accepted,candidates[canI+5,]))
    accepted<-rbind(accepted,candidates[canI+5,]);
  if((!numCan<(canI+6))&&checkFunc(accepted,candidates[canI+6,]))
    accepted<-rbind(accepted,candidates[canI+6,]);
  return(accepted);
}
 PP.fzCaMnRc <- function(accepted,candidates,checkFunc){
  ## First, check alignment, in order to simplify tail call
  numCan <- dim(candidates)[[1]]
  if(numCan<4)
    return(PP.fzCline.addOne(accepted,
                               candidate=candidates,
                               checkFunc=checkFunc));
 
  ## Assumes all intelligent work has been done, so brute force.
  ## Unroll loops, too.
  for(canI in 4*(1:(numCan%/%4))-3){
    if(checkFunc(accepted,candidates[canI,]))
      accepted<-rbind(accepted,candidates[canI,]);
    if(checkFunc(accepted,candidates[canI+1,]))
      accepted<-rbind(accepted,candidates[canI+1,]);
    if(checkFunc(accepted,candidates[canI+2,]))
      accepted<-rbind(accepted,candidates[canI+2,]);
    if(checkFunc(accepted,candidates[canI+3,]))
      accepted<-rbind(accepted,candidates[canI+3,]);
       
  }
  canI= 4*(numCan%/%4)-3
  if((!numCan<(canI+4))&&checkFunc(accepted,candidates[canI+4,]))
    accepted<-rbind(accepted,candidates[canI+4,]);
  if((!numCan<(canI+5))&&checkFunc(accepted,candidates[canI+5,]))
    accepted<-rbind(accepted,candidates[canI+5,]);
  if((!numCan<(canI+6))&&checkFunc(accepted,candidates[canI+6,]))
    accepted<-rbind(accepted,candidates[canI+6,]);
  return(accepted);
} 


PP.fzCline.addBest <- function(accepted,candidates,checkFunc){
  ## This should be called just once per candidate pool.

  ## Useful expressions...

  addFz1 <- PP.fzCline.addOne;
  addFzN <- PP.fzCline.addMany;
  
  ## returns a list of the reduced set and the remaining candidates,
  ## if any. Currently uses mutation.
  pluckSet <- function(pluckTF){
    pluckI <- which(pluckTF)
    pluckO <- which(!(pluckTF)) 
    if(length(pluckI)==0)               #Anything to pluck?
      return()                          #Nope. go away
    if(length(pluckO)==0){
      ## Nothing will be left
      if(length(pluckI)==1){            #Just one? Pull it.
        accepted<<-rbind(accepted,candidates);
        candidates<<-NULL;
        return();
      }
      ## Otherwise, check for duplicates, etc.
      accepted <<- rbind(accepted,
                         addFz1(candidates[1,],
                                candidates[2:(dim(candidates)[[1]]),],
                                checkFunc=checkFunc));
      candidates<<-NULL;
      return();
      ## return(list(rbind(accepted,candidates),NULL));
    }
    ## Something to pluck, something left.
    if(length(pluckI)==1){             #Just one? Pull it.
      accepted  <<- rbind(accepted,candidates[pluckI,]);
      candidates<<- candidates[pluckO,];
      return();
    }
    ## Otherwise, check for duplicates, etc.
    pkFirst <- pluckI[1];
    pluckI <- pluckI[pluckI!=pkFirst]   #There had better not be dupes here...
    accepted   <<- rbind(accepted,
                         addFz1(candidates[pkFirst,],
                                candidates[pluckI,],
                                checkFunc=checkFunc));
    candidates <<- candidates[pluckO,];
    return();
  }

  ## call this to check a channel.  Assumes mutation.
  doPluck <- function(param){
    ## sanity checks
    if(!is.data.frame(candidates))
      return();
    if(!all(param %in% names(accepted),
            param %in% names(candidates)))
      return();
    if(max(accepted[,param]) < (pVal <- max(candidates[,param]))){
      pluckSet(candidates[,param] == pVal);
      if(!is.data.frame(candidates))
        return();
    }
    if(min(accepted[,param]) > (pVal <- min(candidates[,param]))){
      pluckSet(candidates[,param] == pVal);
      if(!is.data.frame(candidates))
        return();
    }
    return();
  }

  ## Assume that all extreme values of the pools are up for consideration
  for(paramI in names(accepted)){
    doPluck(paramI);
  }
  
  if(!is.data.frame(candidates))
    return(accepted);

  return(addFzN(accepted,candidates,checkFunc=checkFunc));
}


PP.fz.oldCheckFunc <- function(paramNames){
  checkCW <- function(accepted,candidate){  
    ## left side bounds
    if(sum(accepted$center<=candidate$center & accepted$width>=candidate$width)==0)
      return(TRUE);
    if(sum(accepted$center>=candidate$center & accepted$width<=candidate$width)==0)
      return(TRUE);
    ##right side bounds
    if(sum(accepted$center>=candidate$center & accepted$width>=candidate$width)==0)
      return(TRUE);
    if(sum(accepted$center<=candidate$center & accepted$width<=candidate$width)==0)
      return(TRUE);
    return(FALSE);
  }
  checkCWpMM <- function(accepted,candidate){  
    ## left side bounds
    if(sum(accepted$center<=candidate$center & accepted$width>=candidate$width&
           accepted$pMin >= candidate$pMin  & accepted$pMax <= candidate$pMax )==0)
      return(TRUE);
    if(sum(accepted$center>=candidate$center & accepted$width<=candidate$width&
           accepted$pMin <= candidate$pMin  & accepted$pMax >= candidate$pMax)==0)
      return(TRUE);
    ##right side bounds
    if(sum(accepted$center>=candidate$center & accepted$width>=candidate$width&
           accepted$pMin >= candidate$pMin  & accepted$pMax <= candidate$pMax )==0)
      return(TRUE);
    if(sum(accepted$center<=candidate$center & accepted$width<=candidate$width&
           accepted$pMin <= candidate$pMin  & accepted$pMax >= candidate$pMax)==0)
      return(TRUE);
    return(FALSE);
  }
  if(identical(sort(paramNames),c("center","width"))){
    return(checkCW);
  }
  if(identical(sort(paramNames),c("center","pMax","pMin","width"))){
    return(checkCWpMM);
  }

  

  ##default, accept all
  return(function(accepted,candidate) return(TRUE));
}

PP.fz.getCheckFunc <- function(paramNames){
  accLeft  <- quote(acc$center<=can$center)
  accRight <- quote(acc$center>=can$center)
  accWIn   <- quote(acc$width >=can$width )
  accWOut  <- quote(acc$width <=can$width )
  accMMOut <- quote(acc$pMin  <=can$pMin  &
                    acc$pMax  >=can$pMax  )
  accMMIn  <- quote(acc$pMin  >=can$pMin  &
                    acc$pMax  <=can$pMax  )
  accLTIn  <- quote(acc$deltaL<=can$deltaL&
                    acc$tauL  <=can$tauL  )
  accLTOut <- quote(acc$deltaL>=can$deltaL&
                    acc$tauL  >=can$tauL  )
  accRTIn  <- quote(acc$deltaR<=can$deltaR&
                    acc$tauR  <=can$tauR )
  accRTOut <- quote(acc$deltaR>=can$deltaR&
                    acc$tauR  >=can$tauR  )
  accMTIn  <- quote(acc$deltaM<=can$deltaM&
                    acc$tauM  <=can$tauM  )
  accMTOut <- quote(acc$deltaM>=can$deltaM&
                    acc$tauM  >=can$tauM  )
  check<- quote(if(all(!((accScale)&(accLoc)))) return(TRUE));
  
 
  checkFunc<-function(acc,can) return(TRUE);
  if((!("center" %in% paramNames))||(!("width"%in%paramNames)))
    return(checkFunc);
  accOut<-accWOut;
  accIn <- accWIn;
  if(("pMin" %in% paramNames)&&("pMax" %in% paramNames)){
    accOut<-substitute(accOld&accEx,list(accOld=accOut,accEx=accMMOut));
    accIn<-substitute(accOld&accEx,list(accOld=accIn,accEx=accMMIn));
  } else if(("pMin" %in%paramNames)&&(!("pMax" %in%paramNames)) ||
            ("pMax" %in%paramNames)&&(!("pMin" %in%paramNames))){
    return(checkFunc);
  }
  if(("deltaM"%in%paramNames)&&("tauM"%in%paramNames)){
    accOut<-substitute(accOld&accEx,list(accOld=accOut,accEx=accMTOut));
    accIn<-substitute(accOld&accEx,list(accOld=accIn,accEx=accMTIn));
  } else if(("deltaM"%in%paramNames)&&(!("tauM"%in%paramNames)) ||
            ("tauM"%in%paramNames)&&(!("deltaM"%in%paramNames))){
    return(checkFunc);
  }
  if(("deltaL"%in%paramNames)&&("tauL"%in%paramNames)){
    accLOut<-substitute(accOld&accEx,list(accOld=accOut,accEx=accLTOut));
    accLIn<-substitute(accOld&accEx,list(accOld=accIn,accEx=accLTIn));
    
  } else if(("deltaL"%in%paramNames)&&(!("tauL"%in%paramNames)) ||
            ("tauL"%in%paramNames)&&(!("deltaL"%in%paramNames))){
    return(checkFunc);
  } else {
    accLOut<-accOut;
    accLIn<-accIn;
  }
  cLO<-substitute(substitute(c,list(accLoc=accRight,accScale=accLOut)),list(c=check));
  cLI<-substitute(substitute(c,list(accLoc=accLeft,accScale=accLIn)),list(c=check));
    
  
  if(("deltaR"%in%paramNames)&&("tauR"%in%paramNames)){
    accROut<-substitute(accOld&accEx,list(accOld=accOut,accEx=accRTOut));
    accRIn<-substitute(accOld&accEx,list(accOld=accIn,accEx=accRTIn));
  } else if(("deltaR"%in%paramNames)&&(!("tauR"%in%paramNames)) ||
            ("tauR"%in%paramNames)&&(!("deltaR"%in%paramNames))){
    return(checkFunc);
  } else {
    accROut<-accOut;
    accRIn<-accIn;
  }
  
    cRO<-substitute(substitute(c,list(accLoc=accLeft,accScale=accROut)),list(c=check));
    cRI<-substitute(substitute(c,list(accLoc=accRight,accScale=accRIn)),list(c=check));
    
  checkE<-do.call(expression,
                  list(eval(cLI),
                       eval(cLO),
                       eval(cRI),
                       eval(cRO),
                       quote(return(FALSE))));
  body(checkFunc)<-as.call(c(as.name("{"),checkE));
  ##print(checkFunc);
  return(checkFunc);
}

  
  

hzar.getLLCutParam <- function(dataGroups,params,cutValue=2){
  params<-as.character(params);
  if(inherits(dataGroups, c("hzar.fitRequest","hzar.dataGroup"))){
    dataGroups<-hzar.fit2DataGroup(dataGroups);
  } else if( is.list(dataGroups)){
    if(is.character(names(dataGroups))){
      return(do.call(rbind,sapply(names(dataGroups),
                                  function(x) hzar.getLLCutParam(dataGroups[[x]],params,cutValue),simplify = FALSE)));
    }
    return(do.call(rbind,lapply(dataGroups,hzar.getLLCutParam,params,cutValue)));
  } else {
    stop("hzar.getLLCutParam does not understand class of dataGroups.");
  }
  
  data.param=dataGroups$data.param[dataGroups$data.LL$model.LL >
    max(dataGroups$data.LL$model.LL -cutValue),];
  tempFunc <- function(x){
    res<-list(min(data.param[[x]]),
              max(data.param[[x]]));
    names(res) <- paste(x,cutValue,"LL",c("Low","High"),sep="");
    return(res);
  }
  return(do.call(data.frame,do.call(c,lapply(params,tempFunc))));
}

hzar.getCredParamRed <- function(dataGroup){
  junk<-hzar.getCredParam(dataGroup);
  junk1<-PP.fzCline.reduce(junk);
  junk1<-PP.fzCline.reduce(junk1);
  res<-hzar.make.fzCline(lapply(1:(dim(junk1)[[1]]),
         function(x) hzar.gen.cline(junk1[x,],
                                    dataGroup)));
  return(res);
}

hzar.make.fzCline <- function(clineList){
  list(clines=clineList,
       listFuncInt=function(xVal,funcList=res$clines){
    
         yList <- as.numeric(lapply(funcList,function(x,u) x$clineFunc(u),u=xVal));
         return(data.frame(x=xVal,yMin=min(yList),yMax=max(yList)));
       },
       fzCline=function(xVals,listFunc=res$listFuncInt){
         xVals <- as.numeric(xVals);
         if(length(xVals>1))
           return(do.call(rbind,lapply(xVals,listFunc)));
         if(length(xVals<1))
           return(numeric(0));
         return(listFunc(xVals));
       })->res;
  return(res);
}

