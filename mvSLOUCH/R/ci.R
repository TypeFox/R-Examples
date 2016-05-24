.calcCI<-function(EvolModel,params,data,designToEstim,lPrecalculates,KnownParams,vVars=NULL,conditional=FALSE,minLogLik=-Inf,conf.level=0.95,t=1,regressCovar=NULL){
    n<-ncol(lPrecalculates$mTreeDist)
    
    ## prepare params for this
    vEstimedPoint<-.ci.vectorizeParams(params$paramPoint,KnownParams)
    
    gradLogLik<-grad(.ci.loglikfunc,vEstimedPoint,method="simple",KnownParams=KnownParams,params=params,data=data,lPrecalculates=lPrecalculates,EvolModel=EvolModel,vVars=vVars,conditional=conditional,minLogLik=minLogLik)
    gradLogLik<-matrix(gradLogLik,ncol=1)
    ## method simple not implemented for hessian Richardson takes too long ... so we need to call Jacobian on grad with simple
    hessLogLik<-jacobian(function(x,KnownParams,params,data,lPrecalculates,EvolModel,vVars,conditional,minLogLik){
	    grad(.ci.loglikfunc,x,method="simple",KnownParams=KnownParams,params=params,data=data,lPrecalculates=lPrecalculates,EvolModel=EvolModel,vVars=vVars,conditional=conditional,minLogLik=minLogLik)
	},vEstimedPoint,method="simple",KnownParams=KnownParams,params=params,data=data,lPrecalculates=lPrecalculates,EvolModel=EvolModel,vVars=vVars,conditional=conditional,minLogLik=minLogLik
    )
    if (length(which(is.nan(hessLogLik)))>0){
	hessLogLik[which(is.nan(hessLogLik))]<-0
	print("WARNING: NaNs in the Hessian changed to 0!!!!")
    }
    if (length(which(is.na(hessLogLik)))>0){
	hessLogLik[which(is.na(hessLogLik))]<-0
	print("WARNING: NAs in the Hessian changed to 0!!!!")
    }
    if (length(which(is.infinite(hessLogLik)))>0){
	hessLogLik[which(is.infinite(hessLogLik))]<-0
	print("WARNING: Infs in the Hessian changed to 0!!!!")
    }

    mA<-gradLogLik%*%t(gradLogLik)
    if (min(abs(diag(qr.R(qr(hessLogLik)))))<1e-07){print("WARNING : Log likelihood surface is very flat, confidence intervals are EXTREMELY APPROXIMATE.");mB<-pseudoinverse(hessLogLik)}
    else{mB<-solve(hessLogLik)}
    vCIs<-qnorm(1-(1-conf.level)/2)*sqrt(diag(mB%*%mA%*%t(mB))) 

    lEstimatedPoint<-.ci.listParams(vEstimedPoint,params$paramPoint,KnownParams)        
    ci.lower<-.ci.listParams(vEstimedPoint-vCIs,params$paramPoint,KnownParams)
    ci.upper<-.ci.listParams(vEstimedPoint+vCIs,params$paramPoint,KnownParams)
    lCI<-vector("list",length(lEstimatedPoint))

    lPoint.lower<-params
    lPoint.upper<-params
    vnames<-c()
    for (i in 1:length(lEstimatedPoint)){ 
	lCI[[i]]<-vector("list",3)
	names(lCI[[i]])<-c("Lower.end","Estimated.Point","Upper.end")
	lCI[[i]]$Lower.end<-ci.lower[[i]]
	lCI[[i]]$Estimated.Point<-lEstimatedPoint[[i]]
	lCI[[i]]$Upper.end<-ci.upper[[i]]
	vnames[i]<-paste(names(lEstimatedPoint)[i],".confidence.interval",sep="")
	lPoint.lower[[which(names(lPoint.lower)==names(lEstimatedPoint)[i])]]<-ci.lower[[i]]
	lPoint.upper[[which(names(lPoint.upper)==names(lEstimatedPoint)[i])]]<-ci.upper[[i]]
    }
    names(lCI)<-vnames

    lCI$lower.summary<-NA
    lCI$lower.summary<-tryCatch({.Params.summary(lPoint.lower,EvolModel,designToEstim,NULL,t,-Inf,0,0,NA,lPrecalculates=list(tree.height=lPrecalculates$tree.height))},error=function(e){print(paste("Error in lower confidence interval calculation",e))})
    lCI$upper.summary<-NA
    lCI$upper.summary<-tryCatch({.Params.summary(lPoint.upper,EvolModel,designToEstim,NULL,t,-Inf,0,0,NA,lPrecalculates=list(tree.height=lPrecalculates$tree.height))},error=function(e){print(paste("Error in lower confidence interval calculation",e))})

## Do eigen CIs ------------------------------------------------------------------
## in code we use HL/hl instead of eig as originally this built confidence intervals for half-lives but it turned out that
## it was stabler and easier to build confidence intervals for the eigenvalues
    if (is.element("A",names(params$paramPoint))){
	eigA<-eigen(params$paramPoint$A)
	P<-eigA$vectors
	vEstimedHLs<-Re(eigA$values)
	gradLogLikHLs<-grad(.ci.HL.loglikfunc,vEstimedHLs,method="simple",P=P,EigA=eigA$values,params=params,data=data,lPrecalculates=lPrecalculates,EvolModel=EvolModel,vVars=vVars,conditional=conditional,minLogLik=minLogLik)

	gradLogLikHLs<-matrix(gradLogLikHLs,ncol=1)
	## method simple not implemented for hessian Richardson takes too long ... so we need to call Jacobian on grad with simple
	hessLogLikHLs<-jacobian(function(x,P,params,data,lPrecalculates,EvolModel,vVars,conditional,minLogLik){
		grad(.ci.HL.loglikfunc,x,method="simple",P=P,EigA=eigA$values,params=params,data=data,lPrecalculates=lPrecalculates,EvolModel=EvolModel,vVars=vVars,conditional=conditional,minLogLik=minLogLik)
	    },vEstimedHLs,method="simple",P=P,params=params,data=data,lPrecalculates=lPrecalculates,EvolModel=EvolModel,vVars=vVars,conditional=conditional,minLogLik=minLogLik
	)
	mAHLs<-gradLogLikHLs%*%t(gradLogLikHLs)
	if (min(abs(diag(qr.R(qr(hessLogLikHLs)))))<1e-07){print("WARNING : Log likelihood surface is very flat, confidence intervals for eigenvalues are EXTREMELY APPROXIMATE.");mBHLs<-pseudoinverse(hessLogLikHLs)}
	else{mBHLs<-solve(hessLogLikHLs)}
	vCIsHLs<-qnorm(1-(1-conf.level)/2)*sqrt(diag(mBHLs%*%mAHLs%*%t(mBHLs))) 
    
	ci.HLs.lower<-vEstimedHLs-vCIsHLs
	ci.HLs.upper<-vEstimedHLs+vCIsHLs
	
	lCI$eigenvalues.confidence.interval$eigenvalues.confidence.interval<-cbind(ci.HLs.lower,vEstimedHLs,ci.HLs.upper)
        colnames(lCI$eigenvalues.confidence.interval$eigenvalues.confidence.interval)<-c("Lower.end","Estimated.eigenvalue","Upper.end")

	lPoint.HL.lower<-params
	lPoint.HL.upper<-params
    
	lCI$eigenvalues$lower.summary<-NA
	lCI$eigenvalues$lower.summary<-tryCatch({.Params.summary(lPoint.HL.lower,EvolModel,designToEstim,NULL,t,-Inf,0,0,NA,lPrecalculates=list(tree.height=lPrecalculates$tree.height))},error=function(e){print(paste("Error in lower confidence interval calculation",e))})
	lCI$eigenvalues$upper.summary<-NA
	lCI$eigenvalues$upper.summary<-tryCatch({.Params.summary(lPoint.HL.upper,EvolModel,designToEstim,NULL,t,-Inf,0,0,NA,lPrecalculates=list(tree.height=lPrecalculates$tree.height))},error=function(e){print(paste("Error in lower confidence interval calculation",e))})
    
	hlMinMax<-matrix(0,nrow=length(vEstimedHLs),ncol=2)
	hlMinMax[,1]<-vEstimedHLs-designToEstim$sigmaRule*vEstimedHLs
	hlMinMax[,2]<-vEstimedHLs+designToEstim$sigmaRule*vEstimedHLs
	hlstep<-2*designToEstim$sigmaRule*vEstimedHLs/10000^(1/length(vEstimedHLs))
	hlgrid<-.generategrid(hlMinMax,hlstep)    
	for (i in 1:ncol(hlgrid)){vSmallHL<-which(abs(hlgrid[,i])<1e-10);if (length(vSmallHL>0)){hlgrid<-hlgrid[-vSmallHL,]}} ## remove small half lives
	hlgrid<-cbind(hlgrid,apply(hlgrid,1,.ci.HL.loglikfunc,P=P,EigA=eigA$values,params=params,data=data,lPrecalculates=lPrecalculates,EvolModel=EvolModel,vVars=vVars,conditional=conditional,minLogLik=minLogLik))
	hlnames<-c()
	for (i in 1:length(vEstimedHLs)){hlnames<-c(hlnames,paste("Eig.dir.",i,sep=""))}
	hlnames<-c(hlnames,"loglik")
	colnames(hlgrid)<-hlnames
	hlgrid<-as.data.frame(hlgrid)
	lCI$eigenvalues$eigenvalues.support.grid<-hlgrid
    }
## -------------------------------------------------------------------------------
    if ((nrow(regressCovar)==0)||(ncol(regressCovar)==0)){regressCovar<-NULL}
    if (!is.null(regressCovar)){
	lCI$regression.summary<-list()
	vRegCIs<-qnorm(1-(1-conf.level)/2)*sqrt(diag(regressCovar))
	currVar<-length(vRegCIs) ## we go from the end as at the beginning we could have some fixed effects
	if (is.element("BX0",names(designToEstim)) && designToEstim$BX0){
	    trueBX0<-params$B%*%params$vX0
	    lCI$regression.summary$BX0.regression.confidence.interval<-matrix(NA,ncol=3,nrow=length(trueBX0))
	    colnames(lCI$regression.summary$BX0.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
	    lCI$regression.summary$BX0.regression.confidence.interval[,"Estimated.Point"]<-trueBX0
	    lCI$regression.summary$BX0.regression.confidence.interval[,"Lower.end"]<-trueBX0-vRegCIs[(currVar-length(trueBX0)+1):currVar]
	    lCI$regression.summary$BX0.regression.confidence.interval[,"Upper.end"]<-trueBX0+vRegCIs[(currVar-length(trueBX0)+1):currVar]
	    currVar<-currVar-length(trueBX0)
	}
	if (is.element("X0",names(designToEstim)) && designToEstim$X0){
	    lCI$regression.summary$X0.regression.confidence.interval<-matrix(NA,ncol=3,nrow=length(params$vX0))
	    colnames(lCI$regression.summary$X0.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
	    lCI$regression.summary$X0.regression.confidence.interval[,"Estimated.Point"]<-params$vX0
	    lCI$regression.summary$X0.regression.confidence.interval[,"Lower.end"]<-params$vX0-vRegCIs[(currVar-length(params$vX0)+1):currVar]
	    lCI$regression.summary$X0.regression.confidence.interval[,"Upper.end"]<-params$vX0+vRegCIs[(currVar-length(params$vX0)+1):currVar]
	    currVar<-currVar-length(params$vX0)
	}
	if (is.element("B",names(designToEstim)) && designToEstim$B){
	    if (ncol(params$B)==1){
	    	lCI$regression.summary$B.regression.confidence.interval<-matrix(NA,ncol=3,nrow=nrow(params$B))
		colnames(lCI$regression.summary$B.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		lCI$regression.summary$B.regression.confidence.interval[,"Estimated.Point"]<-params$B[,1]
		lCI$regression.summary$B.regression.confidence.interval[,"Lower.end"]<-params$B[,1]-vRegCIs[(currVar-nrow(params$B)+1):currVar]
		lCI$regression.summary$B.regression.confidence.interval[,"Upper.end"]<-params$B[,1]+vRegCIs[(currVar-nrow(params$B)+1):currVar]
	    }else{
	    	lCI$regression.summary$B.regression.confidence.interval<-vector("list",3)
		names(lCI$regression.summary$B.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		lCI$regression.summary$B.regression.confidence.interval$Estimated.Point<-params$B
		lCI$regression.summary$B.regression.confidence.interval$Lower.end<-params$B-matrix(vRegCIs[(currVar-nrow(params$B)+1):currVar],nrow=nrow(params$B),ncol=ncol(params$B),byrow=TRUE)
		lCI$regression.summary$B.regression.confidence.interval$Upper.end<-params$B+matrix(vRegCIs[(currVar-nrow(params$B)+1):currVar],nrow=nrow(params$B),ncol=ncol(params$B),byrow=TRUE)
	    }
	    currVar<-currVar-length(c(params$B))
	}
	if (is.element("psi",names(designToEstim)) && designToEstim$psi){
	    if (ncol(params$mPsi)==1){
	    	lCI$regression.summary$mPsi.regression.confidence.interval<-matrix(NA,ncol=3,nrow=nrow(params$mPsi))
		colnames(lCI$regression.summary$mPsi.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		lCI$regression.summary$mPsi.regression.confidence.interval[,"Estimated.Point"]<-params$mPsi[,1]
		lCI$regression.summary$mPsi.regression.confidence.interval[,"Lower.end"]<-params$mPsi[,1]-vRegCIs[(currVar-nrow(params$mPsi)+1):currVar]
		lCI$regression.summary$mPsi.regression.confidence.interval[,"Upper.end"]<-params$mPsi[,1]+vRegCIs[(currVar-nrow(params$mPsi)+1):currVar]
	    }else{
	    	lCI$regression.summary$mPsi.regression.confidence.interval<-vector("list",3)
		names(lCI$regression.summary$mPsi.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		lCI$regression.summary$mPsi.regression.confidence.interval$Estimated.Point<-params$mPsi
		lCI$regression.summary$mPsi.regression.confidence.interval$Lower.end<-params$mPsi-matrix(vRegCIs[(currVar-nrow(params$mPsi)+1):currVar],nrow=nrow(params$mPsi),ncol=ncol(params$mPsi),byrow=FALSE)
		lCI$regression.summary$mPsi.regression.confidence.interval$Upper.end<-params$mPsi+matrix(vRegCIs[(currVar-nrow(params$mPsi)+1):currVar],nrow=nrow(params$mPsi),ncol=ncol(params$mPsi),byrow=FALSE)
	    }
	    currVar<-currVar-length(c(params$mPsi))
	}
	if (is.element("psi0",names(designToEstim)) && designToEstim$psi0){
	    lCI$regression.summary$mPsi0.regression.confidence.interval<-matrix(NA,ncol=3,nrow=nrow(params$mPsi0))
	    colnames(lCI$regression.summary$mPsi0.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
	    lCI$regression.summary$mPsi0.regression.confidence.interval[,"Estimated.Point"]<-params$mPsi0
	    lCI$regression.summary$mPsi0.regression.confidence.interval[,"Lower.end"]<-params$mPsi0-vRegCIs[(currVar-nrow(params$mPsi0)+1):currVar]
	    lCI$regression.summary$mPsi0.regression.confidence.interval[,"Upper.end"]<-params$mPsi0+vRegCIs[(currVar-nrow(params$mPsi0)+1):currVar]
	    currVar<-currVar-length(params$mPsi0)
	}
	if (is.element("Y0",names(designToEstim)) && designToEstim$Y0 && !designToEstim$y0AncState ){
	    lCI$regression.summary$Y0.regression.confidence.interval<-matrix(NA,ncol=3,nrow=length(params$vY0))
	    colnames(lCI$regression.summary$Y0.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
	    lCI$regression.summary$Y0.regression.confidence.interval[,"Estimated.Point"]<-params$vY0
	    lCI$regression.summary$Y0.regression.confidence.interval[,"Lower.end"]<-params$vY0-vRegCIs[(currVar-length(params$vY0)+1):currVar]
	    lCI$regression.summary$Y0.regression.confidence.interval[,"Upper.end"]<-params$vY0+vRegCIs[(currVar-length(params$vY0)+1):currVar]
	    currVar<-currVar-length(params$vY0)
	}
    }
    lCI
}

.ci.vectorizeParams<-function(params,KnownParams){
    vpars<-c()
    for (i in 1:length(params)){if (!is.element(names(params)[i],names(KnownParams))){vpars<-c(vpars,params[[i]])}}
    vpars
}

.ci.listParams<-function(vparams,params,KnownParams){
    lPoint<-list()
    vnames<-c()
    vparIndex<-1
    iel<-1
    for (i in 1:length(params)){
	if (!is.element(names(params)[i],names(KnownParams))){
	    parLen<-length(c(params[[i]]))
	    if (!is.na(params[[i]][1])){lPoint[[iel]]<-matrix(vparams[vparIndex:(vparIndex+parLen-1)],ncol=ncol(params[[i]]),nrow=nrow(params[[i]]))}
	    else{lPoint[[iel]]<-NA}
	    vnames[iel]<-names(params)[i]
	    vparIndex<-vparIndex+parLen
	    iel<-iel+1
	}
    }
    names(lPoint)<-vnames
    lPoint
}

.ci.loglikfunc<-function(vpar,KnownParams,params,data,lPrecalculates,EvolModel,vVars,conditional,minLogLik){
    lEvalParams<-params
    vparIndex<-1
    for (i in 1:length(params$paramPoint)){
	j<-which(names(lEvalParams)==names(params$paramPoint)[i])
	if (!is.element(names(params$paramPoint)[i],names(KnownParams))){
	    parLen<-length(c(params$paramPoint[[i]]))
	    if (!is.na(params$paramPoint[[i]][1])){lEvalParams[[j]]<-matrix(vpar[vparIndex:(vparIndex+parLen-1)],ncol=ncol(params$paramPoint[[i]]),nrow=nrow(params$paramPoint[[i]]))}
	    vparIndex<-vparIndex+parLen
	}else{lEvalParams[[j]]<-KnownParams[[which(names(KnownParams)==names(params$paramPoint)[i])]]}
    } 
    if (is.element("A",names(lEvalParams))){
	lEvalParams$precalcMatrices<-NULL
	lEvalParams$precalcMatrices<-.decompEigenA.S(lEvalParams,lPrecalculates,NULL,toCalc=list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    }
    .calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=lEvalParams,vVars=vVars,conditional=conditional,FALSE,minLogLik=minLogLik)
}

.ci.HL.loglikfunc<-function(vHLs,P,EigA,params,data,lPrecalculates,EvolModel,vVars,conditional,minLogLik){
    if (is.complex(EigA)){vHLs<-complex(real=vHLs,imaginary=Im(EigA))}
    params$A<-Re(P%*%diag(vHLs,ncol(P))%*%solve(P))
    params$precalcMatrices<-NULL
    params$precalcMatrices<-.decompEigenA.S(params,lPrecalculates,NULL,toCalc=list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    .calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=params,vVars=vVars,conditional=conditional,FALSE,minLogLik=minLogLik)
}
