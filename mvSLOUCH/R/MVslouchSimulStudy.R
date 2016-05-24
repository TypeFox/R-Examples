.SimulStudy<-function(PhylTree,modelParams,regimes=NULL,regimes.times=NULL,regimes.types=NULL,n=NULL,maxNicheBranch=0,bShortInternal=FALSE,params=NULL,process="unif",procparams=NULL,M.error=NULL){
    if (is.null(params)){params<-list()}
    if (!is.element("EvolModel",names(params))){params$EvolModel<-"mvslouch"}       
    if(params$EvolModel=="slouch"){params$EvolModel<-"mvslouch";print("Call to slouch not implemented yet, using mvslouch. You might need to run again with correct input structure.")}

    if (is.null(PhylTree)){## no tree need to generate it
	apeRandomTree<-rtree(n)
	if (bShortInternal){
	     vinternalEdge<-which(!is.element(apeRandomTree$edge[,2],1:n))
	     apeRandomTree$edge.length[-vinternalEdge]<-apeRandomTree$edge.length[-vinternalEdge]+100*max(apeRandomTree$edge.length)	    
	}
	PhylTree<-ape2ouch(apeRandomTree)
	PhylTree@nodelabels[1:(PhylTree@nnodes-PhylTree@nterm)]<-1:(PhylTree@nnodes-PhylTree@nterm)	
    }else{n<-PhylTree@nterm}
    
    PhylTree@nodelabels[which(PhylTree@nodelabels=="")]<-1:length(which(PhylTree@nodelabels==""))
    if (is.null(regimes.times)){
    	if (maxNicheBranch>0){regimes.times<-.drawRandomNiches(PhylTree,maxNicheBranch,process,procparams)}
	else{regimes.times<-sapply(PhylTree@epochs,rev,simplify=FALSE)}
    }
    
    root.regime<-NULL
    if (!is.list(regimes)){## The regimes are given as a vector in the ouch tree format --- Change them to a list
	root.regime<-regimes[1]
        vregimes<-regimes
        regimes<-sapply(PhylTree@lineages[PhylTree@term],function(epch,vregimes){
    	    epch<-rev(epch)
    	    as.character(sapply(epch[-1],function(reg,vregimes){vregimes[reg]},vregimes=vregimes,simplify=TRUE))
        },vregimes=vregimes,simplify=TRUE)
    }else{root.regime<-regimes[[1]][1]}
    
    if (is.null(regimes)){regimes<-.drawRandomRegimes(PhylTree,regimes.types,regimes.times,process,procparams)}
    bOKregimes<-TRUE
    if ((length(regimes)!=length(regimes.times))||(length(regimes)!=length(PhylTree@epochs))){
	bOKregimes<-FALSE;print("Wrong number of regimes, regimes.times vectors")
    }
    else{
	vOKtimes<-sapply(1:length(regimes),function(i,regimes,epochs,regimes.times){
    	    bret<-TRUE
            if (length(regimes[[i]])!=(length(regimes.times[[i]])-1)){bret<-FALSE;print(paste("Problem with ",i,"th regimes or regimes.times",sep=""))}
            if (length(regimes.times[[i]])<length(epochs[[i]])){bret<-FALSE;print(paste(i,"th regimes.times does not agree with phylogenetic tree",sep=""))}
            if (length(regimes[[i]])<(length(epochs[[i]])-1)){bret<-FALSE;print(paste(i,"th regimes does not agree with phylogenetic tree",sep=""))}
            bret            
        },regimes=regimes,epochs=PhylTree@epochs,regimes.times=regimes.times)
        if (length(vOKtimes[vOKtimes])!=length(vOKtimes)){bOKregimes<-FALSE}
    }

    if (bOKregimes){
	regimes.types.orig<-c()
        for (i in 1: length(regimes)){regimes.types.orig<-c(regimes.types.orig,unique(regimes[[i]]))}
        regimes.types.orig<-sort(unique(regimes.types.orig))
        regimes.types<-1:length(regimes.types.orig)
        regimes<-sapply(regimes,function(vregs,regimes.types.orig){
    	    sapply(vregs,function(orgreg,regimes.types.orig){which(regimes.types.orig==orgreg)},regimes.types.orig=regimes.types.orig)
        },regimes.types.orig=regimes.types.orig,simplify=FALSE)
    
	 kY<-NA
	 if (params$EvolModel=="bm"){kY<-0;kX<-nrow(modelParams$Sxx)}
         if (params$EvolModel=="slouch"){kY<-1;kX<-nrow(modelParams$Sxx)}
         if (params$EvolModel=="ouch"){kY<-nrow(modelParams$A);kX<-0}
         if (params$EvolModel=="mvslouch"){kY<-nrow(modelParams$A);kX<-nrow(modelParams$Sxx)}         
                     
         if (!is.element("method",names(params))){params$method<-"glsgc"}
         if (params$EvolModel=="bm"){params$method<-"maxlik"}
         if (!is.element("bShouldPrint",names(params))){params$bShouldPrint<-TRUE}
         if (!is.element("tol",names(params))){params$tol<-.set.tol(params$method)}
         if (!is.element("maxIter",names(params))){params$maxIter<-.set.maxiter(params$method)}
         if (!is.element("maxTries",names(params))){params$maxTries<-10}
         if (!is.element("minLogLik",names(params))){params$minLogLik<- -Inf}

         params$EstimParams<-.set.estimparams(params,kY,kX,length(regimes.types))
         if ((params$EvolModel!="bm")&&(params$EstimParams$designToEstim$y0AncState)){
            if (!is.null(root.regime)){params$EstimParams$designToEstim$y0Regime<-which(regimes.types.orig==root.regime)}
             else{
                if (is.element("y0Regime",names(params$EstimParams$designToEstim))){params$EstimParams$designToEstim$y0Regime<-which(regimes.types.orig==params$EstimParams$designToEstim$y0Regime)}       
                else{params$EstimParams$designToEstim$y0Regime<-regimes[[1]][1]} ## need to choose something anyway ...    
            }
    	    modelParams$vY0<-matrix(modelParams$mPsi[,params$EstimParams$designToEstim$y0Regime],ncol=1,nrow=nrow(modelParams$A))
    	    if (!is.null(rownames(modelParams$A))){rownames(modelParams$vY0)<-rownames(modelParams$A)}
	    if ((params$EvolModel=="mvslouch")&&(!params$EstimParams$designToEstim$y0OnlyFixed)){modelParams$vY0<-modelParams$vY0-solve(modelParams$A)%*%modelParams$B%*%modelParams$vX0}
        }
        if (!is.element("TerminalLabels",params$EstimParams)){params$EstimParams$TerminalLabels<-PhylTree@nodelabels[PhylTree@term]}    
    
	x0<-NA
	if((params$EvolModel=="mvslouch")||(params$EvolModel=="slouch")){x0<-rep(0,params$EstimParams$kY+params$EstimParams$kX);x0[1:params$EstimParams$kY]<-modelParams$vY0}
	if(params$EvolModel=="ouch"){x0<-rep(0,params$EstimParams$kY);x0[1:params$EstimParams$kY]<-modelParams$vY0}
	if(params$EvolModel=="bm"){x0<-rep(0,params$EstimParams$kX); x0[1:params$EstimParams$kX]<-modelParams$vX0}
	if(params$EvolModel=="mvslouch"){x0[(params$EstimParams$kY+1):(params$EstimParams$kY+params$EstimParams$kX)]<-modelParams$vX0}

	mSimulatedData<-.simulVasicekProcPhylTree(phyltree=PhylTree,EvolModel=params$EvolModel,modelParams=modelParams,EstimationParams=params,regimes=regimes,regimes.times=regimes.times,Simulparams=NULL,dropInternal=TRUE)$ExtantSample
	if (!is.null(M.error)){
	    V<-matrix(0,ncol=ncol(M.error)*nrow(M.error),nrow=ncol(M.error)*nrow(M.error))
	    diag(V)<-c(t(M.error))
	    U<-rmvnorm(1,mean=rep(0,ncol(M.error)*nrow(M.error)),sigma=V)
	    mSimulatedData<-mSimulatedData+ matrix(U,ncol=ncol(M.error),nrow=nrow(M.error),byrow=TRUE)
	}
	dfSimulatedData<-as.data.frame(mSimulatedData[PhylTree@nodelabels[PhylTree@term],])
	colnames(dfSimulatedData)<-colnames(mSimulatedData)
	mSimulatedData<-cbind(mSimulatedData,"index"=1:nrow(mSimulatedData))
	rownames(dfSimulatedData)<-rownames(mSimulatedData)[mSimulatedData[PhylTree@nodelabels[PhylTree@term],"index"]]
## -----------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
	print("##########################################################################################")
	print("Summary of simulation point")
	if ((params$EvolModel=="mvslouch")||(params$EvolModel=="ouch")||(params$EvolModel=="slouch")){modelParams$method<-"glsgc"}
	if (params$EvolModel=="bm"){params$method<-"maxlik"}
	modelParams2<-modelParams
	modelParams$regimeTimes<-modelParams$regimes.times
	modelParams$regimeTypes<-modelParams$regimes.types
	lPrecalculates<-.calculate.Tree.dists(PhylTree,modelParams$EstimParams$TerminalLabels)
	vSimData<-c(t(dfSimulatedData))
	if (params$EvolModel=="mvslouch"){
	    mX<-dfSimulatedData[,(params$EstimParams$kY+1):(params$EstimParams$kY+params$EstimParams$kX)]
	    if (params$EstimParams$kX==1){mX<-matrix(mX,ncol=1,byrow=TRUE)}                
	    mX<-t(mX)                
	    mXmX0<-mX
	    if (params$EstimParams$designToEstim$UseX0){mXmX0<-mXmX0- matrix(modelParams$vX0,ncol=ncol(mX),nrow=length(modelParams$vX0),byrow=FALSE) }
	    params$EstimParams$Data<-vector("list",2)
	    names(params$EstimParams$Data)<-c("mX","mXmX0")
	    params$EstimParams$Data$mX<-mX
	    params$EstimParams$Data$mXmX0<-mXmX0
	}else{params$EstimParams$Data<-list(mX=NULL,mXmX0=NULL)}
	if ((params$EvolModel=="mvslouch")||(params$EvolModel=="slouch")){modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,params$EstimParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),params$EstimParams$Data$mXmX0)}
	if (params$EvolModel=="ouch"){modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,params$EstimParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)}	
	if (params$EvolModel=="bm"){modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,params$EstimParams$designToEstim,list(bCalcA=FALSE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)}	
	TruePointLogLik<-.calc.phyl.LogLik.traits(vSimData,lPrecalculates=lPrecalculates,params$EvolModel,modelParams=modelParams,vVars=params$EstimParams$vVars,conditional=params$EstimParams$conditional,FALSE,minLogLik=params$minLogLik)     
	TruePointRSS<-.calc.phyl.LogLik.traits(vSimData,lPrecalculates=lPrecalculates,params$EvolModel,modelParams=modelParams,vVars=params$EstimParams$vVars,conditional=params$EstimParams$conditional,TRUE,minLogLik=params$minLogLik)     
	
	npar0<-0
	if (params$EvolModel=="slouch"){npar0<-2}
	if ((params$EvolModel=="slouch")||(params$EvolModel=="ouch")){
	    if (!is.element("A",names(params$EstimParams$KnownParams))){
		if (params$EstimParams$Atype=="SingleValueDiagonal"){npar0<-npar0+1}
		if (params$EstimParams$Atype=="Diagonal"){npar0<-npar0+nrow(modelParams$A)}
		if (params$EstimParams$Atype=="SymmetricPositiveDefinite"){npar0<-npar0+nrow(modelParams$A)*(nrow(modelParams$A)+1)/2}
		if (params$EstimParams$Atype=="Symmetric"){npar0<-npar0+nrow(modelParams$A)*(nrow(modelParams$A)+1)/2}
		if (params$EstimParams$Atype=="TwoByTwo"){npar0<-npar0+4}
		if (params$EstimParams$Atype=="UpperTri"){npar0<-npar0+nrow(modelParams$A)*(nrow(modelParams$A)+1)/2}
		if (params$EstimParams$Atype=="LowerTri"){npar0<-npar0+nrow(modelParams$A)*(nrow(modelParams$A)+1)/2}
		if (params$EstimParams$Atype=="DecomposablePositive"){npar0<-npar0+nrow(modelParams$A)*nrow(modelParams$A)}
		if (params$EstimParams$Atype=="DecomposableNegative"){npar0<-npar0+nrow(modelParams$A)*nrow(modelParams$A)}
		if (params$EstimParams$Atype=="DecomposableReal"){npar0<-npar0+nrow(modelParams$A)*nrow(modelParams$A)}
		if (params$EstimParams$Atype=="Invertible"){npar0<-npar0+nrow(modelParams$A)*nrow(modelParams$A)}
	    }
	    if (!is.element("Syy",names(params$EstimParams$KnownParams))){
		if (params$EstimParams$Syytype=="SingleValueDiagonal"){npar0<-npar0+1}
		if (params$EstimParams$Syytype=="Diagonal"){npar0<-npar0+nrow(modelParams$Syy)}
		if (params$EstimParams$Syytype=="Symmetric"){npar0<-npar0+nrow(modelParams$Syy)*(nrow(modelParams$Syy)+1)/2}
		if (params$EstimParams$Syytype=="UpperTri"){npar0<-npar0+nrow(modelParams$Syy)*(nrow(modelParams$Syy)+1)/2}
		if (params$EstimParams$Syytype=="LowerTri"){npar0<-npar0+nrow(modelParams$Syy)*(nrow(modelParams$Syy)+1)/2}
		if (params$EstimParams$Syytype=="Any"){npar0<-npar0+nrow(modelParams$Syy)*nrow(modelParams$Syy)}
	    }
	}

	lTruePointSummary<-.Params.summary(modelParams,params$EvolModel,params$EstimParams$designToEstim,data=vSimData,t=1,LogLik=TruePointLogLik,n=ncol(lPrecalculates$mSpecDist),npar0=npar0,RSS=TruePointRSS,lPrecalculates=list(tree.height=lPrecalculates$tree.height)) ## no need to build CIs here, true point
	modelParams<-.cleanUpModelParams(modelParams)
	modelParams$TruePointSummary<-lTruePointSummary
	modelParams<-.cleanUpModelParams(modelParams)
	print(modelParams)    
	print("##########################################################################################")
	modelParams<-modelParams2
	rm(modelParams2)
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
	if ((params$EvolModel=="mvslouch")||(params$EvolModel=="slouch")||(params$EvolModel=="ouch")){
	    if (params$EvolModel=="ouch"){params$method<-"gridgls"}
	    if ((params$EvolModel=="mvslouch")||(params$EvolModel=="slouch")){params$method<-"gridigls"}
	    gridPoint<- .par.inv.transform(modelParams,params$EstimParams)
	    params$EstimParams$mGrid<-matrix(gridPoint,nrow=1,ncol=length(gridPoint))
	    colnames(params$EstimParams$mGrid)<-names(gridPoint)
	    TruePointlEstimResults<-.PhyloSDEestim(PhylTree,dfSimulatedData,params$EstimParams$kY,regimes,regimes.times,which(regimes.types.orig==root.regime),NULL,params,M.error)
	    print("##########################################################################################")
	    if (params$EvolModel=="ouch"){print("GLS estimation under true simulation point")}
	    if ((params$EvolModel=="mvslouch")||(params$EvolModel=="slouch")){print("iGLS estimation under true simulation point")}	    
	    print(TruePointlEstimResults)
	    print("##########################################################################################")
	}
	if (params$EvolModel=="bm"){TruePointlEstimResults<-"BM model does not have a regression setup at the moment"}
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------    
	if ((params$EvolModel=="mvslouch")||(params$EvolModel=="slouch")||(params$EvolModel=="ouch")){params$method<-"glsgc"}
	if (params$EvolModel=="bm"){params$method<-"maxlik"}
	print("##########################################################################################")
	lEstimResults<-.PhyloSDEestim(PhylTree,dfSimulatedData,params$EstimParams$kY,regimes,regimes.times,which(regimes.types.orig==root.regime),NULL,params,M.error)
	print("Maximum likelihood found point")
	print(lEstimResults)
	print("##########################################################################################")
#---------------------------------------------------------------------------------------------
	list(lEstimResults,TruePointlEstimResults)
    }else{print("The regimes, regimes.times and PhylTree variables are not consistant. Cannot perform estimation. Please see manual on their format.")}
    
}

.drawRandomRegimes<-function(PhylTree,regimes.types,regimes.times,process="unif",procparams=NULL){
    PhylTree@epochs<-sapply(PhylTree@epochs,rev,simplify=FALSE)
    PhylTree@lineages<-sapply(PhylTree@lineages,rev,simplify=FALSE)
    regimes<-sapply(regimes.times,function(regTim){rep(NA,length(regTim)-1)},simplify=FALSE)
    nNAregs<-sum(sapply(regimes,length,simplify=TRUE))    
    iEmptyRegs<-1    
    while(nNAregs>0){
	vCurrEmpty<-which(is.na(regimes[[iEmptyRegs]]))
	if (length(vCurrEmpty>0)){
	    ## We assume a tree structure here -> WE HAVE TO OPERATE ON BRANCHES
	    if (process=="unif"){regimes[[iEmptyRegs]][vCurrEmpty]<-sample(regimes.types,length(vCurrEmpty),replace=TRUE)}
	    if (process=="Poisson"){
		for (j in sort(vCurrEmpty)){
		    if (j==1){regimes[[iEmptyRegs]][j]<-sample(regimes.types,1,prob=procparams$P[procparams$anc,])}
		    else{regimes[[iEmptyRegs]][j]<-sample(regimes.types,1,prob=procparams$P[regimes[[iEmptyRegs]][j-1],])}
		}
	    }
	    nNAregs<-nNAregs-length(vCurrEmpty)
	    vTimes<-regimes.times[[iEmptyRegs]][union(vCurrEmpty,vCurrEmpty+1)]
	    vAncsAll<-PhylTree@lineages[[iEmptyRegs+PhylTree@nnodes-PhylTree@nterm]]
	    vAncsTimes<-PhylTree@times[PhylTree@lineages[[iEmptyRegs+PhylTree@nnodes-PhylTree@nterm]]]
	    vAncs<-vAncsAll[intersect(which(vAncsTimes>=min(vTimes)),which(vAncsTimes<=max(vTimes)))]
	    if ((iEmptyRegs<PhylTree@nterm)&&(nNAregs>0)){
		for (j in (iEmptyRegs+1):PhylTree@nterm){
		    vCurrEmptyJ<-which(is.na(regimes[[j]]))
			if (length(vCurrEmptyJ>0)){
			    vTimesJ<-regimes.times[[j]][union(vCurrEmptyJ,vCurrEmptyJ+1)]
			    vAncsAllJ<-PhylTree@lineages[[j+PhylTree@nnodes-PhylTree@nterm]]
			    vAncsTimesJ<-PhylTree@times[PhylTree@lineages[[j+PhylTree@nnodes-PhylTree@nterm]]]
			    vAncsJ<-vAncsAllJ[intersect(which(vAncsTimesJ>=min(vTimesJ)),which(vAncsTimesJ<=max(vTimesJ)))]
			    vCommonEmpty<-intersect(vAncs,vAncsJ)
			    if (length(vCommonEmpty)>1){## there has to be a common brach
				for (k in 2:length(vCommonEmpty)){
				    zi<-vCommonEmpty[k-1]
				    zj<-vCommonEmpty[k]
				    vNicheIndexes<-intersect(which(vTimes>=PhylTree@times[zi]),which(vTimes<=PhylTree@times[zj]))
				    vNicheIndexes<-vNicheIndexes[-length(vNicheIndexes)]
				    vNicheIndexesJ<-intersect(which(vTimesJ>=PhylTree@times[zi]),which(vTimesJ<=PhylTree@times[zj]))
				    vNicheIndexesJ<-vNicheIndexesJ[-length(vNicheIndexesJ)]
				    regimes[[j]][vCurrEmptyJ[vNicheIndexesJ]]<-regimes[[iEmptyRegs]][vCurrEmpty[vNicheIndexes]]
				    nNAregs<-nNAregs-length(vNicheIndexesJ)
				}
			    }
			}
		}
	    }    
	}
	iEmptyRegs<-iEmptyRegs+1
    }
    regimes
}

.drawRandomNiches<-function(PhylTree,maxNicheBranch,process="unif",procparams=NULL){
    if (process=="unif"){
	lDivisions<-sapply(PhylTree@ancestors,function(i,maxNicheBranch){
					    rev(sort(runif(sample(1:maxNicheBranch,1))))
		},maxNicheBranch=maxNicheBranch,simplify=FALSE)
    }
    if (process=="Poisson"){
	lDivisions<-sapply(PhylTree@ancestors,function(i,maxNicheBranch,tau){
					    rev(unique(sapply(cumsum(rexp(sample(1:maxNicheBranch,1),tau)),function(x){min(x,1)})))
		},maxNicheBranch=maxNicheBranch,tau=procparams$tau,simplify=FALSE)    
    }
    lDivisions[[1]]<-NA
    regimes.times<-vector("list",PhylTree@term)
    for (i in 1:PhylTree@nterm){
	lineageI<-PhylTree@lineages[[i+PhylTree@nnodes-PhylTree@nterm]]
	lineageI<-lineageI[-length(lineageI)]
	lRegsTmp<-sapply(1:length(lineageI),function(j,lineageI,lDivisions,epoch){
	    regimeB<-rep(NA,length(lDivisions[[lineageI[j]]])+1)
	    epochTime<-epoch[j]-epoch[j+1]
	    regimeB[1]<-epoch[j]
	    for (r in 2:(length(lDivisions[[lineageI[j]]])+1)){
		regimeB[r]<-lDivisions[[lineageI[j]]][r-1]*epochTime+epoch[j+1]	    
	    }
	    regimeB
	},lineageI=lineageI,lDivisions=lDivisions,epoch=PhylTree@epochs[[i]],simplify=FALSE)
	regimes.times[[i]]<-lRegsTmp[[1]]
    	if (length(lRegsTmp)>1){for (j in 2:length(lRegsTmp)){regimes.times[[i]]<-c(regimes.times[[i]],lRegsTmp[[j]])}}
	regimes.times[[i]]<-rev(regimes.times[[i]])
	regimes.times[[i]]<-c(0,regimes.times[[i]])
	regimes.times[[i]]<-unique(regimes.times[[i]]) ## needed for Markov process if it went beyond the branch
    }    
    regimes.times
}
