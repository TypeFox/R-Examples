.beginEstimationParams<-function(EvolModel,EstimationParams,dfData,PhylTree,Merror=NULL){
    if ((EvolModel=="mvslouch")||(EvolModel=="ouch")){
	if (EvolModel=="mvslouch"){
	    mX<-dfData[,(EstimationParams$kY+1):(EstimationParams$kY+EstimationParams$kX)]
	    if (EstimationParams$kX==1){mX<-matrix(mX,ncol=1,byrow=TRUE)}
	    mXbm<-as.matrix(mX)
	    mX<-t(mX)
	    mXbm<-as.data.frame(rbind(matrix(NA,ncol=ncol(mXbm),nrow=PhylTree@nnodes-PhylTree@nterm),mXbm))	    
	    rownames(mXbm)<-1:nrow(mXbm)
	    vPredIndex<-sapply(1:length(PhylTree@term),function(i,kY,kX){((kY+kX)*(i-1)+kY+1):((kY+kX)*i)},kY=EstimationParams$kY,kX=EstimationParams$kX,simplify=TRUE)
	    bmEstim<-.bm.estim(mXbm,PhylTree,Merror[vPredIndex,vPredIndex])
	    EstimationParams$Fixed$Sxx<-bmEstim$Sxx
	    EstimationParams$Fixed$vX0<-bmEstim$vX0

	    mXmX0<-mX
	    if (EstimationParams$designToEstim$UseX0){mXmX0<-mXmX0- matrix(EstimationParams$Fixed$vX0,ncol=ncol(mX),nrow=length(EstimationParams$Fixed$vX0),byrow=FALSE) }  
	    EstimationParams$Data<-vector("list",2)
	    names(EstimationParams$Data)<-c("mX","mXmX0")
	    EstimationParams$Data$mX<-mX
	    EstimationParams$Data$mXmX0<-mXmX0	    
	}
	if (EvolModel=="ouch"){}
    }
    EstimationParams    
}

.EvaluatePoint<-function(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,bFull,lPoint=NULL,calcLogLik=TRUE,bSummarizePoint=FALSE,minLogLik=-Inf,t=1){
    LogLik<- minLogLik
    if (((EvolModel=="mvslouch")||(EvolModel=="ouch"))&&(!bFull)&&(!is.null(lPoint))){
        modelParams$A<-lPoint$A
        modelParams$Syy<-lPoint$Syy
    }

    if (EvolModel=="mvslouch"){
    	    tmpvY0<-NA
    	    if (!EstimationParams$designToEstim$y0){tmpvY0<-EstimationParams$Fixed$vY0}
    	    vAncPsi<-NA
    	    vX0<-NA
	    mPsi<-NA
	    mPsi0<-NA
    	    if (!EstimationParams$designToEstim$psi){mPsi<-EstimationParams$Fixed$mPsi}
	    if (is.element("mPsi0",names(EstimationParams$Fixed))){mPsi0<-EstimationParams$Fixed$mPsi0}
	    if (EstimationParams$designToEstim$y0AncState){
		if (!EstimationParams$designToEstim$psi){
		    vAncPsi<-EstimationParams$Fixed$mPsi[,EstimationParams$designToEstim$y0Regime]
		    if (!EstimationParams$designToEstim$psi0 && is.element("mPsi0",names(EstimationParams$Fixed))){vAncPsi<-vAncPsi+EstimationParams$Fixed$mPsi0}
		}
		if (!EstimationParams$designToEstim$y0OnlyFixed && !EstimationParams$designToEstim$B && EstimationParams$designToEstim$UseX0){vX0<-EstimationParams$Fixed$vX0;y0OnlyFixed<-FALSE}
	    }
    	    bDzeta<-TRUE
	    bKappa<-TRUE
	    if (!EstimationParams$designToEstim$B){bDzeta<-FALSE;bKappa<-FALSE}
    	    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=bDzeta,lexptcalc=TRUE,kappacalc=bKappa,interceptcalc=TRUE),EstimationParams$Data$mXmX0)
	    if (!is.na(modelParams$precalcMatrices[[6]][[1]][1])){
		modelParams$mCovPhyl<-modelParams$precalcMatrices[[6]][[1]]
		modelParams$invSXX<-modelParams$precalcMatrices[[6]][[2]]
	    }
	    if (!is.element("mCovPhyl",names(modelParams))){modelParams$mCovPhyl<-NA}
	    if (!is.element("invSXX",names(modelParams))){modelParams$invSXX<-NA}
	    if (bFull){
		modelParams<-.slouch.mv.iGLS(lPrecalculates,vY,modelParams,EstimationParams$designToEstim,EstimationParams$Data$mX,UnknownIntercept=FALSE,X0=EstimationParams$Fixed$vX0,tol=tol,maxIter=maxIter,ShouldPrint=bShouldPrint)	
		if (is.element("mCovPhyl",names(modelParams))){modelParams$mCovPhyl<-NA}
	    }
    }
    if ((EvolModel=="ouch")||(EvolModel=="bm")){bDzeta<-FALSE;bKappa<-FALSE}

    if (EvolModel=="ouch"){
	    tmpvY0<-NA

	    if (!EstimationParams$designToEstim$y0){tmpvY0<-EstimationParams$Fixed$vY0}
    	    if (!EstimationParams$designToEstim$psi){mPsi<-EstimationParams$Fixed$mPsi}
    	    if (is.element("mPsi0",names(EstimationParams$Fixed))){mPsi0<-EstimationParams$Fixed$mPsi0}
    	    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=TRUE),NULL)
    	    if (bFull){
		if (!is.element("mCovPhyl",names(modelParams))){modelParams$mCovPhyl<-NA}
		modelParams<-.ouch.GLS(lPrecalculates,vY,modelParams,EstimationParams$designToEstim,tol)
		if (is.element("mCovPhyl",names(modelParams))){modelParams$mCovPhyl<-NA}
	    }
    }
    if ((calcLogLik)||(!is.null(data))){LogLik<-.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,FALSE,minLogLik=minLogLik)}
    else{if(calcLogLik){LogLik<-NA}}
    lPointSummary<-NULL
    if (bSummarizePoint){
        n<-nrow(lPrecalculates$mTreeDist)
	if ((EvolModel=="mvslouch")||(EvolModel=="ouch")){kY<-nrow(modelParams$A)}
        if ((EvolModel=="mvslouch")||(EvolModel=="bm")||(EvolModel=="slouch")){kX<-ncol(modelParams$Sxx)}
        if (EvolModel=="slouch"){kY<-1;}
	vVars2<-NULL
	if ((EvolModel=="mvslouch")||(EvolModel=="ouch")||(EvolModel=="slouch")){vVars2<-c()}
	if ((is.null(EstimationParams$predictors)&&((EvolModel=="mvslouch")||(EvolModel=="ouch")||(EvolModel=="slouch")))){
	    if ((EvolModel=="mvslouch")||(EvolModel=="ouch")){vVars2<-1:kY}
	}else{
	    if (!is.null(EstimationParams$predictors)){
		NumVar<-length(data)/n
		vVars2<-setdiff(1:NumVar,EstimationParams$predictors)
	    }	
	}
	if(!is.null(data)){RSS<-.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,TRUE,minLogLik=minLogLik,vVars2=vVars2)}
	else{RSS<-NA}
	if (EstimationParams$calcCI){
    	    if(bShouldPrint){print("Calculating confidence intervals can take very long time")}
              modelParamsTmp<-modelParams
              modelParamsTmp$B<-NA
              if (EvolModel=="bm"){mCovPhyl<-.calc.phyl.cov(NULL,NULL,lPrecalculates$mAncestorTimes,NULL,EvolModel,modelParams)}
              else{mCovPhyl<-.calc.phyl.cov(lPrecalculates$mTreeDist,lPrecalculates$mSpecDist[nrow(lPrecalculates$mSpecDist),],NULL,lPrecalculates$vSpeciesPairs,EvolModel,modelParams)}
 
              if (EvolModel=="mvslouch"){
		modelParams$precalcMatrices[[4]]$lDzetaKappa<-.decompEigenA.S(modelParamsTmp,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=bDzeta,lexptcalc=TRUE,kappacalc=bKappa,interceptcalc=TRUE),EstimationParams$Data$mXmX0)[[4]]$lDzetaKappa
                lDesign<-.slouch.mv.design.matrix.YcT(lPrecalculates$mSpecDist,lPrecalculates$mTreeDist,lPrecalculates$mAncestorTimes,modelParams,EstimationParams$designToEstim,EstimationParams$Data$mX,modelParams$precalcMatrices[[4]]$lDzeta,modelParams$precalcMatrices[[4]]$lDzetaKappa,modelParams$precalcMatrices[[3]]$lexpmtA,modelParams$precalcMatrices[[3]]$lexptjA,X0=EstimationParams$Fixed$vX0,lPrecalculates$invmAncestorTimes)
                SYY<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY)),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY))]
	        SYX<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY)),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX)))]
    	        SXY<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX))),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY))]
                SXX<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX))),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX)))]
                invSXX<-solve(SXX)
                V<-SYY-SYX%*%invSXX%*%SXY 
              }
              if (EvolModel=="ouch"){
		lDesign<-.ouch.design.matrix(lPrecalculates$mSpecDist,lPrecalculates$mTreeDist,modelParams,EstimationParams$designToEstim,lexpmtA=modelParams$precalcMatrices[[3]]$lexpmtA,lexptjA=modelParams$precalcMatrices[[3]]$lexptjA)
	        V<-mCovPhyl
              }
              if (EvolModel=="bm"){V<-mCovPhyl}
              V[which(abs(V)<1e-15)]<-0
              if ((EvolModel=="ouch")||(EvolModel=="mvslouch")){ 
    		modelParams$regressCovar<-matrix(0,ncol=ncol(lDesign$D),nrow=ncol(lDesign$D))
    		if (ncol(lDesign$D)>0){modelParams$regressCovar<-pseudoinverse(t(lDesign$D)%*%solve(V)%*%lDesign$D)}
    	      }
	      lPointSummary<-.Params.summary(modelParams,EvolModel,EstimationParams$designToEstim,data=data,t=t,LogLik=LogLik,n=ncol(lPrecalculates$mSpecDist),npar0=length(modelParams$vPoint),RSS=RSS,lPrecalculates=lPrecalculates,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik)
        }else{lPointSummary<-.Params.summary(modelParams,EvolModel,EstimationParams$designToEstim,data=data,t=t,LogLik=LogLik,n=ncol(lPrecalculates$mSpecDist),npar0=length(modelParams$vPoint),RSS=RSS,lPrecalculates=list(tree.height=lPrecalculates$tree.height))}
    }
    list(modelParams=modelParams,LogLik=LogLik,PointSummary=lPointSummary)
}

.GridEvaluatePointA<-function(EvolModel,vY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint=FALSE){    
    if (EvolModel=="ouch"){
    	modelParams$precalcMatrices[[2]]<-.decompEigenA.S(modelParams,NULL,NA,list(bCalcA=FALSE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)[[2]]
        if (!is.element("mCovPhyl",names(modelParams))){modelParams$mCovPhyl<-NA}
        modelParams<-.ouch.GLS(lPrecalculates,vY,modelParams,EstimationParams$designToEstim,tol)
    }
    if (EvolModel=="mvslouch"){
	tmpvY0<-NA
	if (!EstimationParams$designToEstim$y0){tmpvY0<-EstimationParams$Fixed$vY0}
    	vAncPsi<-NA
    	vX0<-NA
	mPsi<-NA
    	mPsi0<-NA
    	if (EstimationParams$designToEstim$psi){mPsi<-EstimationParams$Fixed$mPsi}
	if (is.element("mPsi0",names(EstimationParams$Fixed))){mPsi0<-EstimationParams$Fixed$mPsi0}
	if (EstimationParams$designToEstim$y0AncState){
	    if (!EstimationParams$designToEstim$psi){
		vAncPsi<-EstimationParams$Fixed$mPsi[,EstimationParams$designToEstim$y0Regime]
	    	if (!EstimationParams$designToEstim$psi0 && is.element("mPsi0",names(EstimationParams$Fixed))){vAncPsi<-vAncPsi+EstimationParams$Fixed$mPsi0}
	    }	
	    if (!EstimationParams$designToEstim$y0OnlyFixed && !EstimationParams$designToEstim$B && EstimationParams$designToEstim$UseX0){vX0<-EstimationParams$Fixed$vX0;y0OnlyFixed<-FALSE}
	}
	lSmatricesIntercept<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=FALSE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=TRUE),EstimationParams$Data$mXmX0)
    	if (!is.na(lSmatricesIntercept[[6]][[1]][1])){
	    modelParams$mCovPhyl<-lSmatricesIntercept[[6]][[1]]
	    modelParams$invSXX<-lSmatricesIntercept[[6]][[2]]
	}
	if (!is.element("mCovPhyl",names(modelParams))){modelParams$mCovPhyl<-NA}
	if (!is.element("invSXX",names(modelParams))){modelParams$invSXX<-NA}	
        modelParams$precalcMatrices[[2]]<-lSmatricesIntercept[[2]]
        modelParams$precalcMatrices[[5]]<-lSmatricesIntercept[[5]]                            
	modelParams<-.slouch.mv.iGLS(lPrecalculates,vY,modelParams,EstimationParams$designToEstim,EstimationParams$Data$mX,UnknownIntercept=FALSE,X0=EstimationParams$Fixed$vX0,tol=tol,maxIter=maxIter,ShouldPrint=bShouldPrint)
    }                                                                                                                                                                                                         
    modelParams
}

