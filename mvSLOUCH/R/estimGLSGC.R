.glsgc.estim<-function(dfData,data,EvolModel,PhylTree,EstimationParams,modelParams,lPrecalculates=NULL,tol=c(0.0001,0.0001),maxIter=c(50,50),bShouldPrint=FALSE,maxTries=10,minLogLik=-Inf){
    vY<-c(t(dfData)[1:EstimationParams$kY,])
    if (bShouldPrint){print("beginEstim")}

    attempt2<-1
    doneOK<-FALSE
    
    while(!doneOK && (attempt2<=maxTries)){
	tryCatch({
	    ## Some error was generated see what can be removed to help
	    if (attempt2>1){EstimationParams$StartPoint[1:length(EstimationParams$StartPoint)]<-rnorm(length(EstimationParams$StartPoint))}
	    EstimationParams<-.beginEstimationParams(EvolModel,EstimationParams,dfData,PhylTree,modelParams$Merror)
	    vEstim<-EstimationParams$StartPoint
	    vNames<-names(EstimationParams$StartPoint)
        
	    ## Calculate model for starting point
	    attempt<-1
	    bShouldTry<-TRUE
	    vEstim1<-vEstim
	    print(paste("Starting point of heuristic search procedure : "))#,vEstim))
	    print(vEstim)
	    while(bShouldTry&&(attempt<=maxTries)){
		if (bShouldPrint){print("evalPoint")}
		modelParams<-.par.transform(vEstim,EstimationParams,modelParams)
		tmpEvaluatedPoint<-.EvaluatePoint(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
		gc(verbose=bShouldPrint)
		modelParams<-tmpEvaluatedPoint$modelParams
		MaxLogLik<-tmpEvaluatedPoint$LogLik
		if (is.infinite(MaxLogLik)||is.nan(MaxLogLik)||(is.na(MaxLogLik))){
		    MaxLogLik<- minLogLik
		    if (attempt==maxTries-1){vEstim[1:length(vEstim)]<-rnorm(length(vEstim))/10}else{vEstim<-jitter(vEstim)}
		    print("Point generated error, randomly changing it. New start point : ")
		    print(vEstim)
		}
		else{
		    bShouldTry<-FALSE
		    MaxSearchPoint<-modelParams
		    MaxSearchPointParams<-vEstim
		    vEstim1<-vEstim*100
		}
		attempt<-attempt+1
	    }
	    doneOK<- TRUE
	},error=function(e){paste("Restarting from new point, caught error:",e)})
	attempt2<-attempt2+1
    }
    if (attempt2 > maxTries){MaxLogLik<- NA}

    iter<-1
    if (is.infinite(MaxLogLik)||is.nan(MaxLogLik)||(is.na(MaxLogLik))){
	print("Cannot start search from this point errors generated")
	print(vEstim)
	iter<-maxIter[1]+1
	LogLik<- minLogLik
    }

## run HeuristicSearch algorithm
    while((.calc.vec.dist(vEstim,vEstim1)>=tol[1])&&(iter<=maxIter[1])){
	## maximize params
	vEstim1<-vEstim
	LogLik<- minLogLik
	if (EstimationParams$maximMethod=="optim"){## maximize likelihood by optim
	    tryCatch({
		if (bShouldPrint){print("optim")}
		OptimRes<-optim(
		    par=vEstim,
	    	    fn=.MinusPhylLogLikFunc,EvolModel=EvolModel,modelParams=modelParams,EstimationParams=EstimationParams,lPrecalculates=lPrecalculates,data=data,vNames=vNames,minLogLik=minLogLik
		)
		gc(verbose=bShouldPrint)
		LogLik<-(-1)*OptimRes$value
		if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
		vEstim<-OptimRes$par	   
		names(vEstim)<-vNames
		if (LogLik>MaxLogLik){
		    MaxLogLik<-LogLik
		    MaxSearchPoint<-modelParams
		    MaxSearchPointParams<-vEstim
		    MaxSearchPoint<-.EvaluatePoint(EvolModel,data,vY,MaxSearchPoint,lPrecalculates,EstimationParams,NA,NA,NA,FALSE,lPoint=.par.transform(vEstim,EstimationParams),TRUE,FALSE,minLogLik=minLogLik)$modelParams
		}
	    },error=function(e){print(paste("Error in optim",e))})
	
	}
	if (EstimationParams$maximMethod=="nlminb"){## maximize likelihood by nlminb
	    tryCatch({
		nlminbRes<-nlminb(
	    	    start=vEstim,
	    	    objective=.MinusPhylLogLikFunc,EvolModel=EvolModel,modelParams=modelParams,EstimationParams=EstimationParams,lPrecalculates=lPrecalculates,data=data,vNames=vNames,minLogLik=minLogLik
		)
		gc(verbose=bShouldPrint)
		LogLik<-(-1)*nlminbRes$objective
		if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
		vEstim<-nlminbRes$par	    
		names(vEstim)<-vNames
		if (LogLik>MaxLogLik){
		    MaxLogLik<-LogLik
		    MaxSearchPoint<-modelParams
		    MaxSearchPointParams<-vEstim
		    MaxSearchPoint<-.EvaluatePoint(EvolModel,NULL,vY,MaxSearchPoint,lPrecalculates,EstimationParams,NA,NA,NA,FALSE,lPoint=.par.transform(vEstim,EstimationParams),FALSE,FALSE,minLogLik=minLogLik)$modelParams		    
		}
	    },error=function(e){print(paste("Error in nlminb",e))})
	}
	if (bShouldPrint){
	    print(paste(EstimationParams$maximMethod,"found"))
	    print(c(vEstim,"LogLik"=LogLik))
	    print(paste("in iteration",iter,"of search."))
	}
	attempt<-1
	bShouldTry<-TRUE
	vEstimOrg<-vEstim
	while(bShouldTry&&(attempt<=maxTries)){
	    modelParams<-.par.transform(vEstim,EstimationParams,modelParams)
    	    if (bShouldPrint){print("evalPoint")}
 	    tmpEvaluatedPoint<-.EvaluatePoint(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
	    gc(verbose=bShouldPrint)
	    modelParams<-tmpEvaluatedPoint$modelParams    	    
	    LogLik<-tmpEvaluatedPoint$LogLik
	    if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
	    if (is.infinite(LogLik)||is.nan(LogLik)||(is.na(LogLik))){
	        LogLik<- minLogLik
		vEstim<-jitter(vEstim)
		print("Point generated error, randomly changing it. Point : ")
		print(vEstim)
	    }
	    else{
		bShouldTry<-FALSE
		if (LogLik>MaxLogLik){
	    	    MaxLogLik<-LogLik
	    	    MaxSearchPoint<-modelParams
		    MaxSearchPointParams<-vEstim
		}
	    }
	    attempt<-attempt+1
	}
	if (is.infinite(LogLik)||is.nan(LogLik)||(is.na(LogLik))){
	    attempt<-1
	    bShouldTry<-TRUE	
	    while(bShouldTry&&(attempt<=maxTries)){	    
	    	if (attempt==maxTries-1){vEstim[1:length(vEstim)]<-rnorm(length(vEstim))/10}else{vEstim<-jitter(MaxSearchPointParams)}		
		modelParams<-.par.transform(vEstim,EstimationParams,modelParams)
		tmpEvaluatedPoint<-.EvaluatePoint(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
    		gc(verbose=bShouldPrint)
    		modelParams<-tmpEvaluatedPoint$modelParams
		LogLik<-tmpEvaluatedPoint$LogLik
		if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
		if (is.infinite(LogLik)||is.nan(LogLik)||(is.na(LogLik))){
		    print("Point generated error, randomly changing it. Point : ")
		    print(vEstim)	
    		    LogLik<- minLogLik
		}
		else{
		    bShouldTry<-FALSE
		    print("Restarting search from jittered best found estimate.")
		    if (LogLik>MaxLogLik){
	    		MaxLogLik<-LogLik
	    		MaxSearchPoint<-modelParams
			MaxSearchPointParams<-vEstim
		    }
	        }	    
	    attempt<-attempt+1
	    }	    
	}
	if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
	if (is.infinite(LogLik)||is.nan(LogLik)||(is.na(LogLik))){
	    print("Cannot search from this point errors generated")
	    vEstim<-vEstimOrg
	    print(vEstim)	
	    iter<-maxIter[1]+1
    	    LogLik<- minLogLik
	}		
	iter<-iter+1	
    }    
    if (bShouldPrint){
        if (.calc.vec.dist(vEstim,vEstim1)>=tol[1]){print("Heuristic search procedure did not reach convergence at point")}
        else{print("Search algorithm converged")}
	print("Summarizing result, this might take a while.")
    }
    
## Evaluate result of heuristic search algorithm    
    MaxLikHeuristicSearch<-vector("list",2)
    names(MaxLikHeuristicSearch)<-c("FinalFound","MaxLikFound")
    MaxLikHeuristicSearch$MaxLikFound<-"Same as final found"
    HeuristicSearchFinalFind<-vector("list",4)
    names(HeuristicSearchFinalFind)<-c("HeuristicSearchPointFinalFind","ParamsInModel","ParamSummary","LogLik")
    HeuristicSearchFinalFind$HeuristicSearchPointFinalFind<-c(vEstim,LogLik)
    names(HeuristicSearchFinalFind$HeuristicSearchPointFinalFind)<-c(names(EstimationParams$StartPoint),"LogLik")
    if (bShouldPrint){print("Found estimate : ");print(HeuristicSearchFinalFind$HeuristicSearchPointFinalFind)}
    HeuristicSearchFinalFind$ParamsInModel<-.par.transform(HeuristicSearchFinalFind$HeuristicSearchPointFinalFind[-length(HeuristicSearchFinalFind$HeuristicSearchPointFinalFind)],EstimationParams,modelParams)
    tryCatch({	
	tmpEvaluatedPoint<-.EvaluatePoint(EvolModel,data,vY,HeuristicSearchFinalFind$ParamsInModel,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
	HeuristicSearchFinalFind$ParamsInModel<-tmpEvaluatedPoint$modelParams	
	HeuristicSearchFinalFind$LogLik<-tmpEvaluatedPoint$LogLik    	

	vVars2<-NULL
        n<-nrow(lPrecalculates$mTreeDist)                
        if ((EvolModel=="mvslouch")||(EvolModel=="ouch")||(EvolModel=="slouch")){vVars2<-c()}
	if ((is.null(EstimationParams$predictors)&&((EvolModel=="mvslouch")||(EvolModel=="ouch")||(EvolModel=="slouch")))){
	    if ((EvolModel=="mvslouch")||((EvolModel=="ouch"))){vVars2<-1:EstimationParams$kY}
	}else{
	    if (!is.null(EstimationParams$predictors)){
		NumVar<-length(data)/n
	        vVars2<-setdiff(1:NumVar,EstimationParams$predictors)
	    }
        }                                                                        
	RSS<-.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=tmpEvaluatedPoint$modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,TRUE,minLogLik=minLogLik,vVars2=vVars2)
	if (EstimationParams$calcCI){
	    if(bShouldPrint){print("Calculating confidence intervals can take very long time")}
	    HeuristicSearchFinalFind$ParamSummary<-.Params.summary(HeuristicSearchFinalFind$ParamsInModel,EvolModel,EstimationParams$designToEstim,data,1,HeuristicSearchFinalFind$LogLik,nrow(dfData),length(vEstim),RSS,lPrecalculates=lPrecalculates,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik)
	}else{HeuristicSearchFinalFind$ParamSummary<-.Params.summary(HeuristicSearchFinalFind$ParamsInModel,EvolModel,EstimationParams$designToEstim,data,1,HeuristicSearchFinalFind$LogLik,nrow(dfData),length(vEstim),RSS,lPrecalculates=list(tree.height=lPrecalculates$tree.height))}
	HeuristicSearchFinalFind$ParamsInModel<-.cleanUpModelParams(HeuristicSearchFinalFind$ParamsInModel)
    },error=function(e){print(paste("Cannot evaluate final found point",e))})
    MaxLikHeuristicSearch$FinalFound<-HeuristicSearchFinalFind
    if (LogLik<MaxLogLik){
	HeuristicSearchMaxFind<-vector("list",4)
	names(HeuristicSearchMaxFind)<-c("HeuristicSearchPointMaxLik","ParamsInModel","ParamSummary","LogLik")
	HeuristicSearchMaxFind$HeuristicSearchPointMaxLik<-c(MaxSearchPointParams,MaxLogLik)
	names(HeuristicSearchMaxFind$HeuristicSearchPointMaxLik)<-c(names(EstimationParams$StartPoint),"LogLik")
	if (bShouldPrint){print("Maxmimum likelihood found estimate : ");print(HeuristicSearchMaxFind$HeuristicSearchPointMaxLik)}
        HeuristicSearchMaxFind$ParamsInModel<-MaxSearchPoint
	HeuristicSearchMaxFind$LogLik<-MaxLogLik
	RSS<-.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=HeuristicSearchMaxFind$ParamsInModel,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,TRUE,minLogLik=minLogLik,vVars2=vVars2)
	if (EstimationParams$calcCI){
	    if(bShouldPrint){print("Calculating confidence intervals can take very long time")}
    	    HeuristicSearchMaxFind$ParamSummary<-.Params.summary(HeuristicSearchMaxFind$ParamsInModel,EvolModel,EstimationParams$designToEstim,data,1,HeuristicSearchMaxFind$LogLik,nrow(dfData),length(vEstim),RSS,lPrecalculates=lPrecalculates,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik) 
	}else{HeuristicSearchMaxFind$ParamSummary<-.Params.summary(HeuristicSearchMaxFind$ParamsInModel,EvolModel,EstimationParams$designToEstim,data,1,HeuristicSearchMaxFind$LogLik,nrow(dfData),length(vEstim),RSS,lPrecalculates=list(tree.height=lPrecalculates$tree.height))}
	HeuristicSearchMaxFind$ParamsInModel<-.cleanUpModelParams(HeuristicSearchMaxFind$ParamsInModel)
	MaxLikHeuristicSearch$MaxLikFound<-HeuristicSearchMaxFind           
    }                        
    MaxLikHeuristicSearch
}
