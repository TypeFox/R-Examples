.simulVasicekProcPhylTree<-function(phyltree,EvolModel,modelParams,EstimationParams=NULL,regimes=NULL,regimes.times=NULL,Simulparams=NULL,dropInternal=TRUE,bAllTrajectories=FALSE,M.error=NULL){
## tree is assumed to be in ouch format
## we don't need  descendent list as we have the tree in ouch format use @lineages

    if (bAllTrajectories && (is.null(Simulparams) || !is.element("step",names(Simulparams)))){
	step<-min(c(0.001,phyltree@depth/1000))
	if (is.null(Simulparams)){Simulparams<-list(step=step)}
	else{Simulparams$step<-step}
    }
## temporary simulation solution --- should consider a full model with this
    bJumpAtNode<-FALSE
    if (is.element("jump",names(Simulparams))){bJumpAtNode<-TRUE}
## -----------------------------------------------------------------------    

    if (is.null(regimes.times)){regimes.times<-sapply(phyltree@epochs,rev,simplify=FALSE)}
    if (is.null(regimes)){regimes<-sapply(regimes.times,function(reg){rep("reg.1",length(reg)-1)},simplify=FALSE)}
    if (!is.list(regimes)){## The regimes are given as a vector in the ouch tree format --- Change them to a list
        vregimes<-regimes
        regimes<-sapply(phyltree@lineages[phyltree@term],function(epch,vregimes){
            epch<-rev(epch)
            as.character(sapply(epch[-1],function(reg,vregimes){vregimes[reg]},vregimes=vregimes,simplify=TRUE))
	},vregimes=vregimes,simplify=TRUE)
    }                                                                 
    regimes.types.orig<-c()
    for (i in 1: length(regimes)){regimes.types.orig<-c(regimes.types.orig,unique(regimes[[i]]))}
    regimes.types.orig<-sort(unique(regimes.types.orig))
    regimes.types<-1:length(regimes.types.orig)
    regimes<-sapply(regimes,function(vregs,regimes.types.orig){
        sapply(vregs,function(orgreg,regimes.types.orig){which(regimes.types.orig==orgreg)},regimes.types.orig=regimes.types.orig)
    },regimes.types.orig=regimes.types.orig,simplify=FALSE)
                                                        
    modelParams$regimes<-regimes
    modelParams$regimeTimes<-regimes.times

    
     params<-list()
     params$EvolModel<-EvolModel
     if (params$EvolModel=="slouch"){params$EvolModel<-"mvslouch";print("Call to slouch not yet implemented, using mvslouch model. You might need to run again with correct input structure.")}  
     if (!is.element("method",names(params))){params$method<-"glsgc"}
     if (params$EvolModel=="bm"){params$method<-"maxlik"}
     if (!is.element("bShouldPrint",names(params))){params$bShouldPrint<-TRUE}
     if (!is.element("tol",names(params))){params$tol<-.set.tol(params$method)}
     if (!is.element("maxIter",names(params))){params$maxIter<-.set.maxiter(params$method)}
     if (!is.element("maxTries",names(params))){params$maxTries<-10}
     if (!is.element("minLogLik",names(params))){params$minLogLik<- -Inf}                                                           

    
    lPrecalculates<-.calculate.Tree.dists(phyltree)
    if ((EvolModel=="bm")||(EvolModel=="bmStep")){
	x0<-modelParams$vX0;kX<-nrow(modelParams$Sxx)
	EstimationParams<-.set.estimparams(params,0,kX,length(regimes.types))
    }
    if ((EvolModel=="ouch")||(EvolModel=="ouchStep")){
	kY<-nrow(modelParams$A)
	x0<-modelParams$vY0
	EstimationParams<-.set.estimparams(params,kY,0,length(regimes.types))
        modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,NA,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
	##invlexpmtA<-sapply(modelParams$precalcMatrices[[3]]$lexpmtA,solve,simplify=FALSE)
    }
    if ((EvolModel=="mvslouch")||(EvolModel=="mvslouchStep")){
	kY<-nrow(modelParams$A)
	kX<-ncol(modelParams$B)
    	x0<-c(modelParams$vY0,modelParams$vX0)
    	EstimationParams<-.set.estimparams(params,kY,kX,length(regimes.types))
        modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
	##invlexpmtA<-sapply(modelParams$precalcMatrices[[3]]$lexpmtA,solve,simplify=FALSE)
    }
    RootNode<-phyltree@root
    mTreeTraject<-matrix(NA,ncol=length(x0),nrow=length(phyltree@nodes))
    rownames(mTreeTraject)<-phyltree@nodelabels
    colnames(mTreeTraject)<-names(x0)
    if (!is.na(x0[1])){ mTreeTraject[RootNode,]<-x0} ## the state of the root is x0    
    if (bJumpAtNode){vJumpNode<-rep(NA,nrow(mTreeTraject))}
    if (bAllTrajectories){
    	if (is.element(EvolModel,c("bm","bmStep","ouch","ouchStep","mvslouch","mvslouchStep"))){ ## THIS NEEDS TO BE THOUGHT OUT IN CASE OF FUTURE DEVELOPMENT
	    lFullTraject<-vector("list",length(phyltree@ancestors)-1)
	    currbranch<-1
	    if (!substr(EvolModel,start=nchar(EvolModel)-3,stop=nchar(EvolModel))=="Step"){
	    	orgEvolModel<-EvolModel
		EvolModel<-paste(EvolModel,"Step",sep="")
	    }
	}else{
	    print("Cannot yet simulate full trajectory for the chosen model")
	    bAllTrajectories<-FALSE
	}
    }else{lFullTraject<-NA}

    for (Term in phyltree@term){## for each terminal node -- defines a lineage
	vTermLineage<-rev(phyltree@lineages[[Term]]) ## we want to reverse it to start from root
	i<-1
	while(!is.na(mTreeTraject[vTermLineage[i],1])){i<-i+1}
	if(i<length(vTermLineage)+1){## just to check if we haven't run over
	    for (j in i:length(vTermLineage)){
    		if (bAllTrajectories){
    		    lFullTraject[[currbranch]]<-vector("list",3)
    		    names(lFullTraject[[currbranch]])<-c("branch","nodenames","trajectory")
    		    lFullTraject[[currbranch]]$branch<-c(vTermLineage[j-1],vTermLineage[j])
    		    lFullTraject[[currbranch]]$nodenames<-phyltree@nodelabels[c(vTermLineage[j-1],vTermLineage[j])]
    		}
    		Xprev<-mTreeTraject[vTermLineage[j-1],]
		timeDiff<-phyltree@times[vTermLineage[j]]-phyltree@times[vTermLineage[j-1]]

## temporary simulation solution --- should consider a full model with this	
		if (bJumpAtNode){		    
		    bDoJump<-FALSE
		    if (Simulparams$jump$jumptype=="ForBoth"){
                        if (is.na(vJumpNode[vTermLineage[j-1]])){
                            if (runif(1)<Simulparams$jump$jumpprob){bDoJump<-TRUE}
                            else{bDoJump<-FALSE;vJumpNode[vTermLineage[j-1]]<-0}
                            #bDoJump<-TRUE
                        }
                    }
		    if (Simulparams$jump$jumptype=="RandomLineage"){
			if (is.na(vJumpNode[vTermLineage[j-1]])){vJumpNode[vTermLineage[j-1]]<-sample(c(1,2),1)}
			if ((vJumpNode[vTermLineage[j-1]]==1)||(vJumpNode[vTermLineage[j-1]]==3)){bDoJump<-TRUE}
			if (vJumpNode[vTermLineage[j-1]]==2){bDoJump<-FALSE;vJumpNode[vTermLineage[j-1]]<-3}
		    }
		    if (Simulparams$jump$jumptype=="JumpWithProb"){
		    	if (is.na(vJumpNode[vTermLineage[j-1]])||(vJumpNode[vTermLineage[j-1]]==1)){if (runif(1)<Simulparams$jump$jumpprob){bDoJump<-TRUE}}
			if ((!is.na(vJumpNode[vTermLineage[j-1]]))&&(vJumpNode[vTermLineage[j-1]]==2)){bDoJump<-FALSE;vJumpNode[vTermLineage[j-1]]<-0}
		    }
		    if (bDoJump){		    	
		    	## add according to distribution
		    	if (Simulparams$jump$jumpdistrib=="Normal"){Xprev<-Xprev+rmvnorm(1,mean=Simulparams$jump$vMean,sigma=Simulparams$jump$mCov)}
			
			## any post-add trajectory corrections
			if (Simulparams$jump$jumptype=="ForBoth"){mTreeTraject[vTermLineage[j-1],]<-Xprev;vJumpNode[vTermLineage[j-1]]<-0}
			if (Simulparams$jump$jumptype=="RandomLineage"){vJumpNode[vTermLineage[j-1]]<-0}
			if (Simulparams$jump$jumptype=="JumpWithProb"){
			    if (is.na(vJumpNode[vTermLineage[j-1]])){vJumpNode[vTermLineage[j-1]]<-1}
			    else{vJumpNode[vTermLineage[j-1]]<-vJumpNode[vTermLineage[j-1]]+1}
			}
		    }
		}
## -----------------------------------------------------------------------    
		if (EvolModel=="slouch"){}## Empty at the moment
		if (EvolModel=="bm"){
		    vMean<-Xprev[1:kX]
		    mCov<-timeDiff*t(modelParams$Sxx)%*%modelParams$Sxx
		    Xdrawn<-rmvnorm(n=1,mean=vMean,sigma=mCov) 
		    Xdrawn<-Xdrawn[nrow(Xdrawn),]
		}
		if (EvolModel=="bmStep"){
			itermNum<-Term-(phyltree@nnodes-phyltree@nterm)
	    		vX0<-Xprev[1:kX]
	    		tmpmodelparams<-modelParams
	    		tmpmodelparams$vX0<-matrix(vX0,ncol=1,nrow=kX)
			Xdrawn<-.bm.simulate(step=ifelse(is.null(Simulparams),timeDiff,Simulparams$step),duration=timeDiff,modelParams=tmpmodelparams,regimes=NULL,regimes.times=NULL,mCov=NULL)
			if (bAllTrajectories){
			    Xdrawn[,1]<-Xdrawn[,1]+phyltree@times[vTermLineage[j-1]]
			    lFullTraject[[currbranch]]$trajectory<-Xdrawn[1:(nrow(Xdrawn)-1),1:ncol(Xdrawn)]
			}
			Xdrawn<-Xdrawn[nrow(Xdrawn)-1,2:ncol(Xdrawn)]
		}
		if (EvolModel=="ouchStep"){
			itermNum<-Term-(phyltree@nnodes-phyltree@nterm)
			vWhichTimes<-intersect(which(regimes.times[[itermNum]]>=phyltree@times[vTermLineage[j-1]]),which(regimes.times[[itermNum]]<=phyltree@times[vTermLineage[j]]))
			regimesCurrTimes<-regimes.times[[itermNum]][vWhichTimes[1:(length(vWhichTimes))]]		
			regimesCurr<-regimes[[itermNum]][vWhichTimes[1:(length(vWhichTimes)-1)]]		
	    		vY0<-Xprev[1:kY]
	    		tmpmodelparams<-modelParams
	    		tmpmodelparams$vY0<-matrix(vY0,ncol=1,nrow=kY)
	    		regimesCurrTimes<-regimesCurrTimes-phyltree@times[vTermLineage[j-1]]
			Xdrawn<-.ouch.simulate(step=ifelse(is.null(Simulparams),timeDiff,Simulparams$step),duration=timeDiff,modelParams=tmpmodelparams,regimes=regimesCurr,regimes.times=regimesCurrTimes,mCov=NULL)			
			if (bAllTrajectories){
			    Xdrawn[,1]<-Xdrawn[,1]+phyltree@times[vTermLineage[j-1]]
			    lFullTraject[[currbranch]]$trajectory<-Xdrawn[1:(nrow(Xdrawn)-1),1:ncol(Xdrawn)]
			}
			Xdrawn<-Xdrawn[nrow(Xdrawn)-1,2:ncol(Xdrawn)]
		}
		if (EvolModel=="mvslouchStep"){
			itermNum<-Term-(phyltree@nnodes-phyltree@nterm)
			vWhichTimes<-intersect(which(regimes.times[[itermNum]]>=phyltree@times[vTermLineage[j-1]]),which(regimes.times[[itermNum]]<=phyltree@times[vTermLineage[j]]))
			regimesCurrTimes<-regimes.times[[itermNum]][vWhichTimes[1:(length(vWhichTimes))]]		
			regimesCurr<-regimes[[itermNum]][vWhichTimes[1:(length(vWhichTimes)-1)]]		
	    		vY0<-Xprev[1:kY]
	    		vX0<-Xprev[(kY+1):(kY+kX)]
	    		tmpmodelparams<-modelParams
	    		tmpmodelparams$vY0<-matrix(vY0,ncol=1,nrow=kY)
	    		tmpmodelparams$vX0<-matrix(vX0,ncol=1,nrow=kX)
	    		regimesCurrTimes<-regimesCurrTimes-phyltree@times[vTermLineage[j-1]]
			Xdrawn<-.mvslouch.simulate(step=ifelse(is.null(Simulparams),timeDiff,Simulparams$step),duration=timeDiff,modelParams=tmpmodelparams,regimes=regimesCurr,regimes.times=regimesCurrTimes,mCov=NULL)
			if (bAllTrajectories){
			    Xdrawn[,1]<-Xdrawn[,1]+phyltree@times[vTermLineage[j-1]]
			    lFullTraject[[currbranch]]$trajectory<-Xdrawn[1:(nrow(Xdrawn)-1),1:ncol(Xdrawn)]
			}
			Xdrawn<-Xdrawn[nrow(Xdrawn)-1,2:ncol(Xdrawn)]
		}

		if ((EvolModel=="ouch")||(EvolModel=="mvslouch")){
		    itermNum<-Term-(phyltree@nnodes-phyltree@nterm)
		    vWhichTimes<-intersect(which(regimes.times[[itermNum]]>=phyltree@times[vTermLineage[j-1]]),which(regimes.times[[itermNum]]<=phyltree@times[vTermLineage[j]]))
	    	    vY0<-Xprev[1:kY]
	    	    mPsi<-modelParams$mPsi
	    	    mPsi0<-modelParams$mPsi0
#	    	    expmtA<-modelParams$precalcMatrices[[3]]$lexpmtA[[itermNum]]
#	    	    exptjA<-modelParams$precalcMatrices[[3]]$lexptjA[[itermNum]][vWhichTimes]

		    expmtAcurr<-.calc.exptA(t=-timeDiff,modelParams$precalcMatrices[[1]])   ##A=modelParams$A) ## correct the mean value structure we are only moving along a single branch
##		    exptAcorr<-expmtAcurr%*%.calc.exptA(t=-phyltree@times[vTermLineage[j-1]],A=modelParams$A)%*%invlexpmtA[[itermNum]]    ## and not through the whole tree, we want the mean at the branch end
##		    exptjA<-sapply(modelParams$precalcMatrices[[3]]$lexptjA[[itermNum]][vWhichTimes],function(mexptjA,exptAcorr){exptAcorr%*%mexptjA},exptAcorr=exptAcorr,simplify=FALSE)

		    exptjA<-sapply(c(regimes.times[[itermNum]][vWhichTimes]-regimes.times[[itermNum]][vWhichTimes[1]]-timeDiff),function(t,precalc){.calc.exptA(t=t,precalc)},precalc=modelParams$precalcMatrices[[1]],simplify=FALSE)		    
	    	    regimesCurr<-regimes[[itermNum]][vWhichTimes[1:(length(vWhichTimes)-1)]]		
	    	    names(regimesCurr)<-NULL
		    if (EvolModel=="ouch"){
#		        vMean<-.calc.mean.ouch.mv(expmtA,vY0,mPsi,mPsi0,exptjA,regimesCurr)
		        vMean<-.calc.mean.ouch.mv(expmtAcurr,vY0,mPsi,mPsi0,exptjA,regimesCurr)
		        mCov<-.calc.cov.ouch.mv(timeDiff,modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])
	    	    }		
		    if (EvolModel=="mvslouch"){## we do not simulate the trajectory -> just draw from appropriate distribution
			vX0<-Xprev[(kY+1):(length(Xprev))]
			A1B<-modelParams$precalcMatrices[[1]]$A1B		
#			vMean<-.calc.mean.slouch.mv(expmtA,A1B,vY0,vX0,mPsi,mPsi0,exptjA,regimesCurr)
    			vMean<-.calc.mean.slouch.mv(expmtAcurr,A1B,vY0,vX0,mPsi,mPsi0,exptjA,regimesCurr)
			mCov<-.calc.cov.slouch.mv(timeDiff,modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])		
		    }
		    
		    Xdrawn<-rmvnorm(n=1,mean=vMean,sigma=mCov) 
		    Xdrawn<-Xdrawn[nrow(Xdrawn),]
		}
		if(bAllTrajectories){currbranch<-currbranch+1}
		mTreeTraject[vTermLineage[j],]<-Xdrawn
	    }
	}else{print("Something is wrong with tree --- seems to be a network")}	
    }    
    if (dropInternal){mTreeTraject[-c(phyltree@term),]<-NA}
    
    if (EvolModel=="bm"){if (!is.null(colnames(modelParams$Sxx))){colnames(mTreeTraject)<-colnames(modelParams$Sxx)}}
    if (EvolModel=="ouch"){if (!is.null(colnames(modelParams$A))){colnames(mTreeTraject)<-colnames(modelParams$A)}}
    if (EvolModel=="mvslouch"){if ((!is.null(colnames(modelParams$A)))&&(!is.null(colnames(modelParams$Sxx)))){colnames(mTreeTraject)<-c(colnames(modelParams$A),colnames(modelParams$Sxx))}}
    
    if (!is.null(M.error)){
	M.error<-.createCovariancematrix(M.error,phyltree@nterm,ncol(mTreeTraject),NULL,"measurement error")
    
	m.errors<-rmvnorm(1,mean=rep(0,ncol(M.error)),sigma=M.error)
	m.errors<-matrix(m.errors,nrow=phyltree@nterm,ncol=ncol(mTreeTraject),byrow=TRUE)
	mTreeTraject[c(phyltree@term),]<-mTreeTraject[(phyltree@term),]+m.errors
    }
    mTreeTraject<-as.data.frame(mTreeTraject)
    simulReturn<-vector("list",3)
    names(simulReturn)<-c("Tree","ExtantSample","FullTrajectory")
    simulReturn$Tree<-phyltree
    simulReturn$ExtantSample<-mTreeTraject
    simulReturn$FullTrajectory<-lFullTraject
    simulReturn
}


drawPhylProcess<-function(PhylTraitProcess,phyltree=NULL,vColours="black",plotlayout=c(1,1),additionalfigs=FALSE,modelParams=NULL,EvolModel=NULL,xlimits=NULL,ylimits=NULL){
    ## prepare data matrix to plot
    mData<-NA
    if (is.list(PhylTraitProcess)){
	if (is.element("FullTrajectory",names(PhylTraitProcess))){
	    lPhylTraject<-PhylTraitProcess$FullTrajectory
	    mData<-lPhylTraject[[1]]$trajectory
	    if (length(lPhylTraject)>1){for (i in 2:length(lPhylTraject)){mData<-rbind(mData,lPhylTraject[[i]]$trajectory)}}
	}
    }else{if (is.matrix(PhylTraitProcess)){mData<-PhylTraitProcess}}
    
    if (!is.na(mData[1])){    
	ntraits<-ncol(mData)-1
	## check if we have the colors
	if ((length(vColours)==0)||(is.na(vColours[1]))||(is.null(vColours))){vColours<-"black"}
	if (length(vColours)<ntraits){vColours<-rep(vColours,length.out=ntraits)}
    
	## prepare plot area
	if (prod(plotlayout)<ntraits){plotlayout<-c(1,ntraits);print("WARNING : Possibly ugly plot layout, change parameter plotlayout!")}
    
	## plot
	par(mfrow=(plotlayout))
	for (i in 2:(ntraits+1)){
	    minmaxx<-NA
	    if (!is.null(xlimits)){
		if ((is.vector(xlimits))&&(length(xlimits)==2)){minmaxx<-xlimits}
		if ((is.list(xlimits))&&(length(xlimits)==1)&&(is.vector(xlimits[[1]]))&&(length(xlimits[[1]])==2)){minmaxx<-xlimits[[1]]}
		if ((is.list(xlimits))&&(length(xlimits)==ntraits)&&(is.vector(xlimits[[i-1]]))&&(length(xlimits[[i-1]])==2)){minmaxx<-xlimits[[i-1]]}
		if ((is.matrix(xlimits))&&(ncol(xlimits)==2)&&(nrow(xlimits)==1)){minmaxx<-xlimits[1,]}
		if ((is.matrix(xlimits))&&(ncol(xlimits)==2)&&(nrow(xlimits)==ntraits)){minmaxx<-xlimits[i-1,]}		
	    }
	    if (is.na(minmaxx)){
		minmaxx<-c(min(mData[,i]),max(mData[,i]))
	    }	    
	    minmaxy<-NA
	    if (!is.null(ylimits)){
		if ((is.vector(ylimits))&&(length(ylimits)==2)){minmaxy<-ylimits}
		if ((is.list(ylimits))&&(length(ylimits)==1)&&(is.vector(ylimits[[1]]))&&(length(ylimits[[1]])==2)){minmaxy<-ylimits[[1]]}
		if ((is.list(ylimits))&&(length(ylimits)==ntraits)&&(is.vector(ylimits[[i-1]]))&&(length(ylimits[[i-1]])==2)){minmaxy<-ylimits[[i-1]]}
		if ((is.matrix(ylimits))&&(ncol(ylimits)==2)&&(nrow(ylimits)==1)){minmaxy<-ylimits[1,]}
		if ((is.matrix(ylimits))&&(ncol(ylimits)==2)&&(nrow(ylimits)==ntraits)){minmaxy<-ylimits[i-1,]}		
	    }
	    if (is.na(minmaxy)){
	    	minmaxy<-c(min(mData[,1]),max(mData[,1]))
	    }	    

	    plot(mData[,i],mData[,1],col=vColours[i-1],pch=19,cex=0.2,main="",xlab="",ylab="",frame.plot=FALSE,axes=FALSE,ylim=rev(minmaxy),xlim=minmaxx)
	    if (additionalfigs && !is.null(modelParams) && !is.null(EvolModel)){
		if (EvolModel=="bm"){abline(v=modelParams$vX0[i-1,1],lty=2,lwd=1.5);mtext(side=3,at=modelParams$vX0[i-1,1],text=expression(X[0]),cex=2)}
		if (EvolModel=="ouou"){abline(v=modelParams$vY0[i-1,1],lty=2,lwd=1.5);abline(v=modelParams$mPsi[i-1,1],lty=2,lwd=1.5);mtext(side=3,at=modelParams$vY0[i-1,1],text=expression(X[0]),cex=2);mtext(side=3,at=modelParams$mPsi[i-1,1],text=expression(theta),cex=2)}
		if (EvolModel=="mvslouch"){
		    kY<-nrow(modelParams$A);kX<-nrow(modelParams$Sxx)
		    if (i-1<=kY){abline(v=modelParams$vY0[i-1,1],lty=2,lwd=1.5);abline(v=modelParams$mPsi[i-1,1],lty=2,lwd=1.5);mtext(side=3,at=modelParams$vY0[i-1,1],text=expression(Y[0]),cex=2);mtext(side=3,at=modelParams$mPsi[i-1,1],text=expression(psi),cex=2)}
		    else{abline(v=modelParams$vX0[i-1-kY,1],lty=2,lwd=1.5);mtext(side=3,at=modelParams$vX0[i-1-kY,1],text=expression(X[0]),cex=2)}
		}
	    }
	}
    }else{print("Error: wrong data provided to plotting function")}
    NA        
}

