.Params.summary<-function(modelParams,EvolModel,designToEstim,data=NULL,t=1,LogLik=-Inf,n=0,npar0=0,RSS=NA,lPrecalculates=NULL,KnownParams=NULL,conf.level=0.95,vVars=NULL,conditional=FALSE,minLogLik=-Inf){
       tryCatch({
	modelParams$designToEstim<-designToEstim
	tree.height<-NA
	if ((!is.null(lPrecalculates))&&(!is.null(lPrecalculates$tree.height))){tree.height<-lPrecalculates$tree.height}
	names(LogLik)<-c()
	lParamSummary=switch(EvolModel,
    		bm=.params.summary.bm(modelParams,data,LogLik,RSS,n,tree.height),
        	ouch=.params.summary.ouch(modelParams,data,t,LogLik,n,npar0,RSS,tree.height),
#        	slouch=.params.summary.slouch(modelParams,data,LogLik,n,npar0,RSS,tree.height),
        	mvslouch=.params.summary.mvslouch(modelParams,data,t,LogLik,n,npar0,RSS,tree.height)
    	    )
	if ((!is.null(lPrecalculates))&&(!is.null(lPrecalculates$mSpecDist))){
	    regressCovar<-NULL
	    if (is.element("regressCovar",names(modelParams))){regressCovar<-modelParams$regressCovar}
	    modelParams$paramPoint<-.cleanUpModelParams(modelParams) 
    	    tryCatch({
    		lParamSummary$confidence.interval<-.calcCI(EvolModel,modelParams,data,designToEstim,lPrecalculates,KnownParams,vVars,conditional,minLogLik,conf.level,t,regressCovar)
	    },error=function(e){print(paste("Error in confidence interval calculation",e))})
	}
	lParamSummary
    },error=function(e){print(paste("Error in parameter summary",e))})
}

.params.summary.bm<-function(modelParams,data,LogLik,RSS,n,tree.height){
    lParamsSummary<-vector("list",1)
    names(lParamsSummary)<-c("StS")
    lParamsSummary$StS<-modelParams$Sxx%*%t(modelParams$Sxx)
    numobs<- n*ncol(modelParams$Sxx)
    if (!is.null(data)){numobs<- n*ncol(modelParams$Sxx)-length(which(is.na(c(data))))}
    lParamsSummary$LogLik<-LogLik
    lParamsSummary$dof<- ncol(modelParams$Sxx)*(ncol(modelParams$Sxx)+1)/2+ncol(modelParams$Sxx)
    lParamsSummary$m2loglik<- -2*LogLik
    lParamsSummary$aic<- -2*LogLik+2*lParamsSummary$dof
    lParamsSummary$aic.c<- lParamsSummary$aic +2*lParamsSummary$dof*(lParamsSummary$dof+1)/(numobs-lParamsSummary$dof-1)
    lParamsSummary$sic<- lParamsSummary$m2loglik+log(numobs)*lParamsSummary$dof
    lParamsSummary$bic<-lParamsSummary$m2loglik+lParamsSummary$dof*log(numobs)
    lParamsSummary$RSS<-RSS
    lParamsSummary
}

.params.summary.ouch<-function(modelParams,data=NULL,t=1,LogLik=-Inf,n=0,npar0=0,RSS=NA,tree.height=1){
    kY<-ncol(modelParams$A)
    lParamsSummary<-vector("list",18)   
    names(lParamsSummary)<-c("phyl.halflife","expmtA","mPsi.rotated","mPsi0.rotated","cov.matrix","corr.matrix","trait.regression","stationary.cov.matrix","stationary.corr.matrix","StS","LogLik","dof","m2loglik","aic","aic.c","sic","bic","RSS")
    lDecomps<-.decompEigenA.S(modelParams,NULL,NA,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    lParamsSummary$expmtA<-.calc.exptA(-t,lDecomps[[1]])
    lParamsSummary$mPsi.rotated<-apply(modelParams$mPsi,2,function(vPsi,expmtA){(diag(1,nrow(expmtA),ncol(expmtA))-expmtA)%*%vPsi},expmtA=lParamsSummary$expmtA)
    lParamsSummary$mPsi0.rotated<-(diag(1,nrow(lParamsSummary$expmtA),ncol(lParamsSummary$expmtA))-lParamsSummary$expmtA)%*%modelParams$mPsi0
    lParamsSummary$cov.matrix<-.calc.cov.ouch.mv(t,lDecomps[[1]],lDecomps[[2]])
    lParamsSummary$corr.matrix<-cov2cor(lParamsSummary$cov.matrix)
    if(kY>1){
	lParamsSummary$trait.regression<-NULL
	tryCatch({
	    lParamsSummary$trait.regression<-sapply(1:kY,function(i,mCov){mCov[i,-i]%*%solve(mCov[-i,-i])},mCov=lParamsSummary$cov.matrix,simplify=FALSE)
	},error=function(e){print(paste("Error in trait regression calculation",e))})
    }else{lParamsSummary$trait.regression<-NULL}
    lParamsSummary$phyl.halflife<-.calc.phyl.halflife(modelParams$A,tree.height)
    k<-ncol(modelParams$A)
    if (length(which(Re(lDecomps[[1]]$eigA$values)<=0))==0){
	hadInvL1pL2<-apply(matrix(0:(k^2-1),k,k,byrow=TRUE),c(1,2),.CalcVlqStat,vlambda=lDecomps[[1]]$eigA$values,k=k)
	lParamsSummary$stationary.cov.matrix<-Re(lDecomps[[1]]$eigA$vectors%*%
					    (hadInvL1pL2*(
						lDecomps[[1]]$invP%*%(lDecomps[[2]]$S11)%*%t(lDecomps[[1]]$invP)))%*%
					    t(lDecomps[[1]]$eigA$vectors))
	lParamsSummary$stationary.corr.matrix<-cov2cor(lParamsSummary$stationary.cov.matrix)
    }else{
        lParamsSummary$stationary.cov.matrix<-NULL
        lParamsSummary$stationary.cov.matrix.comment<-"A has negative eigenvalues, stationary covariance does not exist"
        lParamsSummary$stationary.corr.matrix<-NULL
        lParamsSummary$stationary.corr.matrix.comment<-"A has negative eigenvalues, stationary correlation does not exist"
    }
    lParamsSummary$StS<-lDecomps[[2]]$S11
    lParamsSummary$LogLik<-LogLik
    lParamsSummary$dof<-npar0
    if (modelParams$designToEstim$psi){lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)*ncol(modelParams$mPsi)}
    if (modelParams$designToEstim$psi0){lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)}
    if (!modelParams$designToEstim$y0AncState && modelParams$designToEstim$y0){lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)}
    lParamsSummary$m2loglik<- -2*LogLik
    numNAdata<-0
    if (!is.null(data)){numNAdata<-length(which(is.na(data)))}
    lParamsSummary$aic<-lParamsSummary$m2loglik+2*lParamsSummary$dof
    lParamsSummary$aic.c<-lParamsSummary$aic+2*lParamsSummary$dof*(lParamsSummary$dof+1)/(nrow(modelParams$A)*n-numNAdata-lParamsSummary$dof-1)
    lParamsSummary$sic<-lParamsSummary$m2loglik+log(nrow(modelParams$A)*n-numNAdata)*lParamsSummary$dof
    lParamsSummary$bic<-lParamsSummary$m2loglik+lParamsSummary$dof*log(nrow(modelParams$A)*n-numNAdata)
    lParamsSummary$RSS<-RSS
    lParamsSummary
}

.params.summary.mvslouch<-function(modelParams,data=NULL,t=1,LogLik=-Inf,n=0,npar0=0,RSS=NA,tree.height=1){
    lParamsSummary<-vector("list",26)   
    names(lParamsSummary)<-c("phyl.halflife","expmtA","optimal.regression","mPsi.rotated","mPsi0.rotated","cov.matrix","corr.matrix","conditional.cov.matrix","conditional.corr.matrix","stationary.cov.matrix","stationary.corr.matrix","optima.cov.matrix","optima.corr.matrix","cov.with.optima","corr.with.optima","evolutionary.regression","trait.regression","StS","LogLik","dof","m2loglik","aic","aic.c","sic","bic","RSS")
    lDecomps<-.decompEigenA.S(modelParams,NULL,NA,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    lParamsSummary$expmtA<-.calc.exptA(-t,lDecomps[[1]])
    lParamsSummary$mPsi.rotated<-apply(modelParams$mPsi,2,function(vPsi,expmtA){(diag(1,nrow(expmtA),ncol(expmtA))-expmtA)%*%vPsi},expmtA=lParamsSummary$expmtA)
    lParamsSummary$mPsi0.rotated<-(diag(1,nrow(lParamsSummary$expmtA),ncol(lParamsSummary$expmtA))-lParamsSummary$expmtA)%*%modelParams$mPsi0
    lParamsSummary$cov.matrix<-.calc.cov.slouch.mv(t,lDecomps[[1]],lDecomps[[2]])
    lParamsSummary$corr.matrix<-cov2cor(lParamsSummary$cov.matrix)
    lParamsSummary$phyl.halflife<-.calc.phyl.halflife(modelParams$A,tree.height)
    kY<-ncol(modelParams$A)
    kX<-ncol(modelParams$B)
    lParamsSummary$evolutionary.regression<-lParamsSummary$cov.matrix[1:kY,(kY+1):(kY+kX)]%*%solve(lParamsSummary$cov.matrix[(kY+1):(kY+kX),(kY+1):(kY+kX)])
    lParamsSummary$trait.regression<-NULL
    tryCatch({
	lParamsSummary$trait.regression<-sapply(1:kY,function(i,mCov){mCov[i,-i]%*%solve(mCov[-i,-i])},mCov=lParamsSummary$cov.matrix,simplify=FALSE)
    },error=function(e){print(paste("Error in trait regression calculation",e))})
    lParamsSummary$conditional.cov.matrix<-lParamsSummary$cov.matrix[1:kY,1:kY]-lParamsSummary$cov.matrix[1:kY,(kY+1):(kY+kX)]%*%solve(lParamsSummary$cov.matrix[(kY+1):(kY+kX),(kY+1):(kY+kX)])%*%lParamsSummary$cov.matrix[(kY+1):(kY+kX),1:kY]
    lParamsSummary$conditional.corr.matrix<-cov2cor(lParamsSummary$conditional.cov.matrix)
    if (length(which(Re(lDecomps[[1]]$eigA$values)<=0))==0){
	hadInvL1pL2<-apply(matrix(0:(kY^2-1),kY,kY,byrow=TRUE),c(1,2),.CalcVlqStat,vlambda=lDecomps[[1]]$eigA$values,k=kY)
	lParamsSummary$stationary.cov.matrix<-Re(lDecomps[[1]]$eigA$vectors%*%
					    (hadInvL1pL2*(
						lDecomps[[1]]$invP%*%(
						    lDecomps[[2]]$S11+
						    lDecomps[[2]]$S12%*%t(lDecomps[[1]]$A1B)+
						    lDecomps[[1]]$A1B%*%lDecomps[[2]]$S21+
				    		    lDecomps[[1]]$A1B%*%lDecomps[[2]]$S22%*%t(lDecomps[[1]]$A1B)
						)%*%t(lDecomps[[1]]$invP)))%*%t(lDecomps[[1]]$eigA$vectors))
        lParamsSummary$stationary.corr.matrix<-cov2cor(lParamsSummary$stationary.cov.matrix)
        lParamsSummary$optimal.regression<- (-1)*lDecomps[[1]]$A1B
    }else{
        lParamsSummary$stationary.cov.matrix<-NULL
        lParamsSummary$stationary.cov.matrix.comment<-"A has negative eigenvalues, stationary covariance does not exist"
        lParamsSummary$stationary.corr.matrix<-NULL
        lParamsSummary$stationary.corr.matrix.comment<-"A has negative eigenvalues, stationary correlation does not exist"
        lParamsSummary$optimal.regression<- NULL
        lParamsSummary$optimal.regression.comment<- "A has negative eigenvalues, optimal regression does not exist"
        lParamsSummary$A1B<-lDecomps[[1]]$A1B
    }
    lParamsSummary$optima.cov.matrix<-t*lDecomps[[1]]$A1B%*%lDecomps[[2]]$S22%*%t(lDecomps[[1]]$A1B)
    lParamsSummary$optima.corr.matrix<-cov2cor(lParamsSummary$optima.cov.matrix)
    lParamsSummary$cov.with.optima<-lParamsSummary$cov.matrix[1:kY,(kY+1):(kY+kX)]%*%t(lDecomps[[1]]$A1B)*(-1) ## minus one here !!
    lParamsSummary$corr.with.optima<-apply(matrix(0:((kY)^2-1),kY,kY,byrow=TRUE),c(1,2),function(ij,kY,mcov.traits,mcov.optima,mcov.with){i<-ij%/%kY+1;j<-ij%%kY+1;mcov.with[i,j]/(sqrt(mcov.traits[i,i]*mcov.optima[j,j]))},kY=kY,mcov.traits=lParamsSummary$cov.matrix,mcov.optima=lParamsSummary$optima.cov.matrix,mcov.with=lParamsSummary$cov.with.optima)
    lParamsSummary$StS<-rbind(cbind(lDecomps[[2]]$S11,lDecomps[[2]]$S12),cbind(lDecomps[[2]]$S21,lDecomps[[2]]$S22))
    lParamsSummary$LogLik<-LogLik
    lParamsSummary$dof<-npar0 + nrow(modelParams$Sxx)+nrow(modelParams$Sxx)*(1+nrow(modelParams$Sxx))/2 ## Sxx is for now estimated directly
    if (modelParams$designToEstim$B){lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$B)*ncol(modelParams$B)}    
    if (modelParams$designToEstim$psi){lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)*ncol(modelParams$mPsi)}    
    if (modelParams$designToEstim$psi0){lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)}    
    if (!modelParams$designToEstim$y0AncState && modelParams$designToEstim$y0){lParamsSummary$dof<- lParamsSummary$dof+nrow(modelParams$A)}
    lParamsSummary$m2loglik<- -2*LogLik
    lParamsSummary$aic<-lParamsSummary$m2loglik+2*lParamsSummary$dof
    numNAdata<-0
    if (!is.null(data)){numNAdata<-length(which(is.na(data)))}
    lParamsSummary$aic.c<-lParamsSummary$aic+2*lParamsSummary$dof*(lParamsSummary$dof+1)/((nrow(modelParams$A)+nrow(modelParams$Sxx))*n-numNAdata-lParamsSummary$dof-1)
    lParamsSummary$sic<-lParamsSummary$m2loglik+log((nrow(modelParams$A)+nrow(modelParams$Sxx))*n-numNAdata)*lParamsSummary$dof
    lParamsSummary$bic<-lParamsSummary$m2loglik+lParamsSummary$dof*log((nrow(modelParams$A)+nrow(modelParams$Sxx))*n-numNAdata)
    lParamsSummary$RSS<-RSS
    lParamsSummary
}

.norm.max<-function(M){max(abs(M))}

.calc.phyl.halflife<-function(A,tree.height){
    mPhylHalfLife<-matrix(NA,nrow=3,ncol=nrow(A))
    rownames(mPhylHalfLife)<-c("eigenvalues","halflife","%treeheight")
    eigA<-eigen(A)
    mPhylHalfLife[1,]<-eigA$values
    mPhylHalfLife[2,]<- log(2)/(Re(eigA$values)) 
    mPhylHalfLife[3,]<- 100*(mPhylHalfLife[2,]/tree.height)
    ## Idea for bound taken from von Lohan, Sensitivity of the matrix exponential remove det part
    list("directions"=eigA$vectors,"halflives"=mPhylHalfLife,"halflifeLowerbounds"=c(Re(log(2)/sum(eigA$values))))
}

.SummarizeFullPoint<-function(vPoint,dfData,PhylTree,EvolModel,regimes.times,regimes,EstimationParams=NULL,modelParams=NULL,t=c(1),dof=NULL,calcCI=NULL,sigmaRule=NULL,tol=0.0001,maxIter=50,bShouldPrint=FALSE,Merror=NULL,predictors=NULL,minLogLik=-Inf,LogLik=NULL){
    if (!is.null(Merror)){
	Merror<-.createCovariancematrix(Merror,nrow(dfData),ncol(dfData),NULL,"measurement error in summary function")
	##if (ncol(Merror)==ncol(dfData)){
	##    tmpMerror<-matrix(0,ncol=ncol(dfData)*nrow(dfData),nrow=ncol(dfData)*nrow(dfData))
	##    diag(tmpMerror)<-c(t(Merror))
	##    Merror<-tmpMerror
	##}
    }
    regimeTimes<-regimes.times
    if (is.null(EstimationParams)){EstimationParams<-list()}
    if (!is.null(calcCI)){
	EstimationParams$calcCI<-calcCI
	if (!is.null(sigmaRule)){
	    if (!is.element("designToEstim",names(EstimationParams))){EstimationParams$desginToEstim<-list()}
	    EstimationParams$desginToEstim$sigmaRule<-sigmaRule
	}
    }
    if (EvolModel=="slouch"){EvolModel<-"mvslouch";print("Call to slouch not yet implemented, using mvslouch model. You might need to run again with correct input structure.")}  
    if((!(is.null(dfData)))&&(!(is.null(PhylTree)))){
	PhylTree@nodelabels[which(PhylTree@nodelabels=="")]<-1:length(which(PhylTree@nodelabels==""))
	if(!is.data.frame(dfData)){
    	    dfData<-as.data.frame(dfData)
    	    rownames(dfData)<-PhylTree@nodelabels[PhylTree@term]
	}
	if (is.null(rownames(dfData))){rownames(dfData)<-PhylTree@nodelabels[PhylTree@term]}    
	colsdfData<-colnames(dfData)
	dfData<-dfData[PhylTree@nodelabels[PhylTree@term],]
    
	## the following is needed if we have one variable previous operation will change dfData to a vector
	if(!is.data.frame(dfData)){ dfData<-as.data.frame(dfData) ;   rownames(dfData)<-PhylTree@nodelabels[PhylTree@term]  ;    colnames(dfData)<-colsdfData }
        if (is.null(regimeTimes)){regimeTimes<-sapply(PhylTree@epochs,rev,simplify=FALSE)}
	if (is.null(regimes)){regimes<-sapply(regimeTimes,function(reg){rep("reg.1",length(reg)-1)},simplify=FALSE)}
	root.regime<-NULL
	if (!is.list(regimes)){## The regimes are given as a vector in the ouch tree format --- Change them to a list
	    root.regime<-regimes[1]
    	    vregimes<-regimes
    	    regimes<-sapply(PhylTree@lineages[PhylTree@term],function(epch,vregimes){
    			    epch<-rev(epch)
    			    as.character(sapply(epch[-1],function(reg,vregimes){vregimes[reg]},vregimes=vregimes,simplify=TRUE))
    	    },vregimes=vregimes,simplify=TRUE)
	}                                                       
	regimes.types.orig<-c()
	for (i in 1: length(regimes)){regimes.types.orig<-c(regimes.types.orig,unique(regimes[[i]]))}
	regimes.types.orig<-sort(unique(regimes.types.orig))
	regimeTypes<-1:length(regimes.types.orig)
	regimes<-sapply(regimes,function(vregs,regimes.types.orig){
    		sapply(vregs,function(orgreg,regimes.types.orig){which(regimes.types.orig==orgreg)},regimes.types.orig=regimes.types.orig)
	    },regimes.types.orig=regimes.types.orig,simplify=FALSE)
    }else{regimeTypes<-c("reg.1");root.regime<-1}

    if (!is.element("kX",names(EstimationParams))){
	    EstimationParams$kX=switch(EvolModel,
    		bm=nrow(modelParams$Sxx),
        	ouch=0,
        	slouch=nrow(modelParams$Sxx),
        	mvslouch=nrow(modelParams$Sxx)
    	    )    
    }
    if (!is.element("kY",names(EstimationParams))){
	EstimationParams$kY=switch(EvolModel,
    		bm=0,
        	ouch=nrow(modelParams$A),
        	slouch=1,
        	mvslouch=nrow(modelParams$A)
    	    )
    }

    EstimationParams<-.set.estimparams(params=list(method="glsgc",EvolModel=EvolModel,EstimParams=EstimationParams,tol=tol,maxIter=maxIter,bShouldPrint=bShouldPrint,Merror=Merror),kY=EstimationParams$kY,kX=EstimationParams$kX,numregs=length(regimeTypes))                                                     
    if (!is.null(predictors)){
	EstimationParams$predictors<-predictors
        predictors<-colnames(dfData)[predictors]
    }
                                            
    if(!is.null(PhylTree)){lPrecalculates<-.calculate.Tree.dists(PhylTree,UserTermLabels=EstimationParams$TerminalLabels)}

    if (EstimationParams$designToEstim$y0AncState){
         if (!is.null(root.regime)){EstimationParams$designToEstim$y0Regime<-which(regimes.types.orig==root.regime)}
         else{
            if (is.element("y0Regime",names(EstimationParams$designToEstim))){EstimationParams$designToEstim$y0Regime<-which(regimes.types.orig==EstimationParams$designToEstim$y0Regime)}       
            else{EstimationParams$designToEstim$y0Regime<-regimes[[1]][1]} ## need to choose something anyway ...    
        }
    }
                                                                       
    
    ## dfData is assumed to be in the following format
    ## number of rows is the number of species, first columns 1:kY response, columns (kY+1):(kY+kX) predictor (for mvslouch)
    ## needs to be done with t as c(M[,cols]) vectorizes M[,cols] by column
    if (!is.null(dfData)){
	data<-c(t(dfData))
	vY<-c(t(dfData)[1:EstimationParams$kY,])            
    }else{data<-NULL;vY<-NA}
    
    if (is.null(vPoint)){
        if (EvolModel=="mvslouch"){
                if(is.element("Fixed",names(EstimationParams))&& is.element("Sxx",names(EstimationParams$Fixed))){Sxx<-EstimationParams$Fixed$Sxx}
                if(is.element("Fixed",names(EstimationParams))&& is.element("vX0",names(EstimationParams$Fixed))){vX0<-EstimationParams$Fixed$vX0}
                if(is.element("Fixed",names(EstimationParams))&& is.element("B",names(EstimationParams$Fixed))){B<-EstimationParams$Fixed$B}
        }
        if ((EvolModel=="ouch")||(EvolModel=="mvslouch")){
            if(is.element("Fixed",names(EstimationParams))&& is.element("vY0",names(EstimationParams$Fixed))){vY0<-EstimationParams$Fixed$vY0}
            if(is.element("Fixed",names(EstimationParams))&& is.element("mPsi",names(EstimationParams$Fixed))){mPsi<-EstimationParams$Fixed$mPsi}
            if(is.element("Fixed",names(EstimationParams))&& is.element("mPsi0",names(EstimationParams$Fixed))){mPsi0<-EstimationParams$Fixed$mPsi0}
        }
    }

    EstimationParams<-.beginEstimationParams(EvolModel,EstimationParams,dfData,PhylTree,Merror)
    if (is.null(vPoint)){
        if (!is.element("Fixed",names(EstimationParams))){EstimationParams$Fixed<-list()}
        if (EvolModel=="mvslouch"){
            if (!is.element("Sxx",names(EstimationParams$Fixed))){EstimationParams$Fixed$Sxx<-Sxx}else{EstimationParams$Fixed$Sxx<-modelParams$Sxx}
	    if (!is.element("vX0",names(EstimationParams$Fixed))){EstimationParams$Fixed$vX0<-vX0}else{EstimationParams$Fixed$vX0<-modelParams$vX0}
    	    if (!is.element("B",names(EstimationParams$Fixed))){EstimationParams$Fixed$B<-B}else{EstimationParams$Fixed$B<-modelParams$B}
        }
        if ((EvolModel=="ouch")||(EvolModel=="mvslouch")){
            if (!is.element("vY0",names(EstimationParams$Fixed))){EstimationParams$Fixed$vY0<-vY0}else{EstimationParams$Fixed$vY0<-modelParams$vY0}
            if (!is.element("mPsi",names(EstimationParams$Fixed))){EstimationParams$Fixed$mPsi<-mPsi}else{EstimationParams$Fixed$mPsi<-modelParams$mPsi}
            if (!is.element("mPsi0",names(EstimationParams$Fixed))){EstimationParams$Fixed$mPsi0<-mPsi0}else{EstimationParams$Fixed$mPsi0<-modelParams$mPsi0}
        }
    }
    if (is.null(modelParams)){
	modelParams<-.par.transform(vPoint,EstimationParams)
	modelParams$vPoint<-vPoint
	bFull<-TRUE
    }else{bFull<-FALSE}
    modelParams$regimeTimes<-regimeTimes
    modelParams$regimes<-regimes
    modelParams$regimeTypes<-regimeTypes
    lEvalPoint<-vector("list",length(t))
    if (is.null(vPoint)){
	if (is.element("A",names(modelParams))){EstimationParams$Fixed$A<-modelParams$A}
	if (is.element("Syy",names(modelParams))){EstimationParams$Fixed$Syy<-modelParams$Syy}
	if (is.null(LogLik)){lEvalPoint<-sapply(t,function(tcalc){.EvaluatePoint(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,FALSE,list(A=EstimationParams$Fixed$A,Syy=EstimationParams$Fixed$Syy),TRUE,TRUE,minLogLik=minLogLik,t=tcalc)},simplify=FALSE)}
	else{lEvalPoint<-sapply(t,function(tcalc){.EvaluatePoint(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,FALSE,list(A=EstimationParams$Fixed$A,Syy=EstimationParams$Fixed$Syy),FALSE,TRUE,minLogLik=LogLik,t=tcalc)},simplify=FALSE)}
    }
    else{
	if (is.null(LogLik)){lEvalPoint<-sapply(t,function(tcalc){.EvaluatePoint(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,bFull,NULL,TRUE,TRUE,minLogLik=minLogLik,t=tcalc)},simplify=FALSE)}
	else{lEvalPoint<-sapply(t,function(tcalc){.EvaluatePoint(EvolModel,data,vY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,bFull,NULL,FALSE,TRUE,minLogLik=LogLik,t=tcalc)},simplify=FALSE)}
    }
    names(lEvalPoint)<-t
    for (i in 1:length(t)){
	lEvalPoint[[i]]$modelParams<-.cleanUpModelParams(lEvalPoint[[i]]$modelParams)
	lEvalPoint[[i]]<-.correct.names(lEvalPoint[[i]],regimes.types.orig,if(EvolModel!="bm"){colnames(dfData)[1:EstimationParams$kY]}else{NULL},if ((EvolModel=="mvslouch")||(EvolModel=="bm")||(EvolModel=="slouch")){colnames(dfData)[(EstimationParams$kY+1):(EstimationParams$kY+EstimationParams$kX)]}else{NULL},predictors,EvolModel)
    }
    if (!is.null(dof)){
	n<-nrow(dfData)
        numNAdata<-length(which(is.na(dfData)))
	for (i in 1:length(t)){
	    lEvalPoint[[i]]$PointSummary$dof<-dof
	    lEvalPoint[[i]]$PointSummary$aic<-lEvalPoint[[i]]$PointSummary$m2loglik+2*dof
	    lEvalPoint[[i]]$PointSummary$aic.c<-lEvalPoint[[i]]$PointSummary$aic+2*dof*(dof+1)/(ncol(dfData)*n-numNAdata-dof-1)
	    lEvalPoint[[i]]$PointSummary$sic<-lEvalPoint[[i]]$PointSummary$m2loglik+log(ncol(dfData)*n-numNAdata)*dof
	    lEvalPoint[[i]]$PointSummary$bic<-lEvalPoint[[i]]$PointSummary$sic
	}
    }
    lEvalPoint
}

.cleanUpModelParams<-function(modelParams){
    if (is.element("eigenSignsA",names(modelParams))){modelParams$eigenSignsA<-NULL}
    if (is.element("GivensQCsignsA",names(modelParams))){modelParams$GivensQCsignsA<-NULL}
    if (is.element("regimes",names(modelParams))){modelParams$regimes<-NULL}
    if (is.element("regimeTypes",names(modelParams))){modelParams$regimeTypes<-NULL}
    if (is.element("regimeTimes",names(modelParams))){modelParams$regimeTimes<-NULL}
    if (is.element("regimes.types",names(modelParams))){modelParams$regimes.types<-NULL}
    if (is.element("regimes.times",names(modelParams))){modelParams$regimes.times<-NULL}
    if (is.element("mCovPhyl",names(modelParams))){modelParams$mCovPhyl<-NULL}
    if (is.element("invSXX",names(modelParams))){modelParams$invSXX<-NULL}
    if (is.element("intercept",names(modelParams))){modelParams$intercept<-NULL}
    if (is.element("precalcMatrices",names(modelParams))){modelParams$precalcMatrices<-NULL}
    if (is.element("paramPoint",names(modelParams))){modelParams$paramPoint<-NULL}    
    if (is.element("regressCovar",names(modelParams))){modelParams$regressCovar<-NULL}    
    if (is.element("EstimParams",names(modelParams))){modelParams$EstimParams<-NULL}    
    if (is.element("kY",names(modelParams))){modelParams$kY<-NULL}    
    if (is.element("kX",names(modelParams))){modelParams$kX<-NULL}    
    if (is.element("method",names(modelParams))){modelParams$method<-NULL}        
    if (is.element("tol",names(modelParams))){modelParams$tol<-NULL}    
    if (is.element("maxIter",names(modelParams))){modelParams$maxIter<-NULL}    
    if (is.element("bShouldPrint",names(modelParams))){modelParams$bShouldPrint<-NULL}    
    if (is.element("EvolModel",names(modelParams))){modelParams$EvolModel<-NULL}    
    if (is.element("maxTries",names(modelParams))){modelParams$maxTries<-NULL}    
    if (is.element("process",names(modelParams))){modelParams$process<-NULL}    
    if (is.element("procparams",names(modelParams))){modelParams$procparams<-NULL}    
    if (is.element("minLogLik",names(modelParams))){modelParams$minLogLik<-NULL}    
    if (is.element("designToEstim",names(modelParams))){modelParams$designToEstim<-NULL}    
    if (is.element("Merror",names(modelParams))){modelParams$Merror<-NULL}    
    modelParams
}
                                                                                                                                                                                                                                                        