.calcLogLikParametrized<-function(vparams,data,vY,PhylTree,EvolModel,lPrecalculates,EstimationParams,modelParams,bFull,tol,maxIter,bShouldPrint,minLogLik=-Inf,parNames=NULL){
    names(vparams)<-parNames
    modelParams<-.par.transform(vparams,EstimationParams,modelParams)
    .EvaluatePoint(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,bFull,NULL,TRUE,FALSE,minLogLik)$LogLik
}

.MinusPhylLogLikFunc<-function(par,EvolModel,modelParams,EstimationParams,lPrecalculates,data,vNames,minLogLik=-Inf){
    LogLik<- minLogLik
    tryCatch({
	    if (EvolModel=="ouch"){
		LogLik<-.MinusPhylLogLikFuncouch(par,modelParams,EstimationParams,lPrecalculates,data,vNames,minLogLik)
	    }
	    if (EvolModel=="mvslouch"){
		LogLik<-.MinusPhylLogLikFuncMVslouch(par,modelParams,EstimationParams,lPrecalculates,data,vNames,minLogLik)
	    }
	},error=function(e){print(paste("Caught:",e))})
    LogLik
}

.MinusPhylLogLikFuncMVslouch<-function(par,modelParamsHeuristicSearch,EstimationParams,lPrecalculates,data,vNames,minLogLik=-Inf){
    names(par)<-vNames
    modelParams<-.par.transform(par,EstimationParams,modelParamsHeuristicSearch)
    modelParams$B<-modelParamsHeuristicSearch$B
    modelParams$mPsi<-modelParamsHeuristicSearch$mPsi
    modelParams$mPsi0<-modelParamsHeuristicSearch$mPsi0
    modelParams$vY0<-modelParamsHeuristicSearch$vY0
    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    (-1)*.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,"mvslouch",modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,FALSE,minLogLik=minLogLik)
}

.MinusPhylLogLikFuncouch<-function(par,modelParamsHeuristicSearch,EstimationParams,lPrecalculates,data,vNames,minLogLik=-Inf){
    names(par)<-vNames
    modelParams<-.par.transform(par,EstimationParams,modelParamsHeuristicSearch)
    modelParams$mPsi<-modelParamsHeuristicSearch$mPsi
    modelParams$mPsi0<-modelParamsHeuristicSearch$mPsi0
    modelParams$vY0<-modelParamsHeuristicSearch$vY0
    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    (-1)*.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,"ouch",modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,FALSE,minLogLik=minLogLik)                
}

.calc.phyl.LogLik.traits<-function(data,lPrecalculates,EvolModel,modelParams,vVars=NULL,conditional=FALSE,RSS=FALSE,minLogLik=-Inf,vVars2=NULL){
    phylLogLik<-minLogLik  
    if (((EvolModel!="ouch")&&(EvolModel!="mvslouch"))||(!(is.na(modelParams$precalcMatrices[[1]]$invA[1])))){
	tryCatch({
	    if(!(is.matrix(lPrecalculates$mSpecDist))){lPrecalculates$mSpecDist<-matrix(lPrecalculates$mSpecDist,nrow=1)}
	    vMean<-.calc.phyl.mean(lPrecalculates$mSpecDist,EvolModel,modelParams)
	    if (is.element("mCovPhyl",names(modelParams))&&!is.na(modelParams$mCovPhyl)){mCov<-modelParams$mCovPhyl}
	    else{mCov<-.calc.phyl.cov(lPrecalculates$mTreeDist,lPrecalculates$mSpecDist[nrow(lPrecalculates$mSpecDist),],lPrecalculates$mAncestorTimes,lPrecalculates$vSpeciesPairs,EvolModel,modelParams)}
	    if (!(is.null(vVars))){
    		n<-ncol(lPrecalculates$mSpecDist)
		kYkX<-length(vMean)/n
		if ((max(vVars)>kYkX)||(min(vVars)<1)){print(paste("Provided variables out bound 1:",kYkX," correcting to this",sep=""))}
		vVars<-intersect(unique(vVars),1:kYkX) ## we are not really interested in degenerate disitributions as same likelihood and singularity in dmvnorm
		vElemsToKeeps<-sort(c(sapply(vVars,function(x){seq(from=x,length.out=n,by=kYkX)})))
		if (length(vElemsToKeeps)<length(vMean)){## no point if we are taking the same variables since we don't consider degenerate ones
		    if (conditional){## conditional on the other variables distribution
			vMean<-vMean[vElemsToKeeps]+mCov[vElemsToKeeps,-vElemsToKeeps]%*%pseudoinverse(mCov[-vElemsToKeeps,-vElemsToKeeps])%*%(data[-vElemsToKeeps]-vMean[-vElemsToKeeps])
	    		mCov<-mCov[vElemsToKeeps,vElemsToKeeps]-mCov[vElemsToKeeps,-vElemsToKeeps]%*%pseudoinverse(mCov[-vElemsToKeeps,-vElemsToKeeps])%*%mCov[-vElemsToKeeps,vElemsToKeeps]
		    }else{## just the marginal distribution of the variables
			vMean<-vMean[vElemsToKeeps]
			mCov<-mCov[vElemsToKeeps,vElemsToKeeps]
		    }
		}
		data<-data[vElemsToKeeps]		
	    }
	    numpoints<-length(vMean)
	    vNAdata<-which(is.na(data))## there are missing values we take the marginal likelihood
	    if (!RSS){
		if (length(vNAdata)>0){data<-data[-vNAdata];vMean<-vMean[-vNAdata];mCov<-mCov[-vNAdata,-vNAdata]}
		phylLogLik<-dmvnorm(x=data, mean=vMean, sigma=mCov, log=TRUE)
		if (phylLogLik<minLogLik){phylLogLik<-minLogLik}
	    }else{
		if (!is.null(vVars2)){
		    phylLogLik<-vector("list",2)
		    names(phylLogLik)<-c("RSS","R2")
		    n<-ncol(lPrecalculates$mSpecDist)
		    kYkX<-numpoints/n
		    vVars2<-intersect(1:kYkX,unique(vVars2)) ## we are not really interested in degenerate disitributions as same likelihood and singularity in dmvnorm
		    kY<-length(vVars2)
		    if (kY<kYkX){
			vElemsToKeeps2<-sort(c(sapply(vVars2,function(x){seq(from=x,length.out=n,by=kYkX)})))
			vElemsToKeeps2mNA<-union(vElemsToKeeps2,vNAdata)
			vElemsToKeeps2<-setdiff(vElemsToKeeps2,vNAdata)			
			vMean2<-vMean[vElemsToKeeps2]+mCov[vElemsToKeeps2,-vElemsToKeeps2mNA]%*%pseudoinverse(mCov[-vElemsToKeeps2mNA,-vElemsToKeeps2mNA])%*%(data[-vElemsToKeeps2mNA]-vMean[-vElemsToKeeps2mNA])
	    		mCov2<-mCov[vElemsToKeeps2,vElemsToKeeps2]-mCov[vElemsToKeeps2,-vElemsToKeeps2mNA]%*%pseudoinverse(mCov[-vElemsToKeeps2mNA,-vElemsToKeeps2mNA])%*%mCov[-vElemsToKeeps2mNA,vElemsToKeeps2]
		    	data2<-data[vElemsToKeeps2]		
		    }else{
			if (length(vNAdata)>0){data<-data[-vNAdata];vMean<-vMean[-vNAdata];mCov<-mCov[-vNAdata,-vNAdata]}
			vMean2<-vMean
	    		mCov2<-mCov
		    	data2<-data		    
		    }
		    mDesign<-matrix(1,ncol=1,nrow=n)%x%diag(1,kY,kY)
		    if (kY<kYkX){
			mDesign<-matrix(1,ncol=1,nrow=n)%x%diag(1,kYkX,kYkX)
			if(ncol(mDesign)>1){mDesign<-mDesign[vElemsToKeeps2,]}else{mDesign<-matrix(c(mDesign)[vElemsToKeeps2],ncol=1)}
		    }else{
			if (length(vNAdata)>0){if(ncol(mDesign)>1){mDesign<-mDesign[-vNAdata,]}else{mDesign<-matrix(c(mDesign)[-vNAdata],ncol=1)}}
		    }
		    vIntercp<- pseudoinverse(t(mDesign)%*%pseudoinverse(mCov2)%*%mDesign)%*%t(mDesign)%*%pseudoinverse(mCov2)%*%data2
		    vMean3<-matrix(1,ncol=1,nrow=n)%x%vIntercp
		    if (kY<kYkX){
			vMean3<-vMean3[vElemsToKeeps2]
		    }else{if (length(vNAdata)>0){vMean3<-vMean3[-vNAdata]}}
		    RSS3<-(t(data2-vMean3))%*%pseudoinverse(mCov2)%*%(data2-vMean3)
		    RSS2<-(t(data2-vMean2))%*%pseudoinverse(mCov2)%*%(data2-vMean2)
	    	    if (length(vNAdata)>0){data<-data[-vNAdata];vMean<-vMean[-vNAdata];mCov<-mCov[-vNAdata,-vNAdata]}
	    	    phylLogLik$RSS<-(data-vMean)%*%pseudoinverse(mCov)%*%(data-vMean)
		    phylLogLik$R2<-1-RSS2/RSS3
		}
		else{phylLogLik<-(data-vMean)%*%pseudoinverse(mCov)%*%(data-vMean)}
	    }
	    },warning=function(w){print(paste("Warning in dmvnorm",w))},error=function(e){print(e)}
	)
    }
    phylLogLik
}


.min.mLogLik.BmPsi<-function(par,X,Y,V){(-1)*dmvnorm(x=Y,mean=X%*%par,sigma=V,log=TRUE)}
    

.min.RSS.BmPsi<-function(par,X,Y,estimParams,Bcols,modelParams,mTreeDist,mSpecDist,vSpeciesPairs,invSXX,kY,kX){
    estimParams<-par
    B<-matrix(estimParams[Bcols],nrow=kY,ncol=kX,byrow=TRUE)
    modelParams$B<-B
    modelParams$precalcMatrices[[1]]$A1B<-modelParams$precalcMatrices[[1]]$invA%*%B
    lCovars<-.slouch.mv.residual.covar(modelParams,B,mTreeDist,mSpecDist[nrow(mSpecDist),],vSpeciesPairs,NA,invSXX)
    V<-lCovars[[1]]
    invV<-solve(V)
    t(Y-X%*%estimParams)%*%invV%*%(Y-X%*%estimParams)
}

                            