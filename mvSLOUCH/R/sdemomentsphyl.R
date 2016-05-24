.calc.phyl.mean<-function(mSpecDist,EvolModel,modelParams){
    vMean=switch(EvolModel,
	bm=.bm.phyl.mean(modelParams$vX0,ncol(mSpecDist)), 
	ouch=.ouch.phyl.mean(mSpecDist,modelParams), 
	slouch=.mvslouch.phyl.mean(mSpecDist,modelParams),
	mvslouch=.mvslouch.phyl.mean(mSpecDist,modelParams)
    )
    vMean[which(abs(vMean)<1e-15)]<-0
    vMean
}

.calc.phyl.cov<-function(mTreeDist,vSpeciesTime,mAncestorTimes,vSpeciesPairs,EvolModel,modelParams){
    mCov=switch(EvolModel,
	bm=.bm.phyl.cov(mAncestorTimes,S=modelParams$Sxx), 
	ouch=.ouch.phyl.cov(mTreeDist,vSpeciesTime,modelParams,vSpeciesPairs), 
	slouch=.mvslouch.phyl.cov(mTreeDist,vSpeciesTime,modelParams,vSpeciesPairs),
	mvslouch=.mvslouch.phyl.cov(mTreeDist,vSpeciesTime,modelParams,vSpeciesPairs)
    )
    if (is.element("Merror",names(modelParams))){mCov<-mCov+modelParams$Merror}
    mCov[which(abs(mCov)<1e-15)]<-0
    mCov
}

.bm.phyl.mean<-function(vY0,n){rep(vY0,n)} 

.bm.phyl.cov<-function(mAncestorTimes,StS=NULL,S=NULL){
    if (is.null(StS)){StS<-S%*%t(S)}
    mCov<-mAncestorTimes%x%StS
    mCov
}

.ouch.phyl.mean<-function(mSpecDist,modelParams){
    if (is.null(modelParams$precalcMatrices[[3]])){
	lexpmtA<-sapply(mSpecDist[nrow(mSpecDist),],function(t){.calc.exptA(-t,modelParams$precalcMatrices[[1]])},simplify=FALSE) ## last row in mSpecDist will be the times of current species
        if (!(is.null(modelParams$regimeTimes))){lexptjA<-
    	    sapply(1:length(modelParams$regimeTimes),function(i,regimeTimes,vSpecDist){tjs<-regimeTimes[[i]];specT<-vSpecDist[i];sapply(tjs,function(t,specT){.calc.exptA(t-specT,modelParams$precalcMatrices[[1]])},specT=specT,simplify=FALSE)},regimeTimes=modelParams$regimeTimes,vSpecDist=mSpecDist[nrow(mSpecDist),],simplify=FALSE)}
        else{ lexptjA<-sapply(1:ncol(mSpecDist),function(i,k){list(lexpmtA[[i]],diag(1,nrow=k,ncol=k))},k=length(modelParams$vY0),simplify=FALSE) }
    }else{
	lexpmtA<-modelParams$precalcMatrices[[3]]$lexpmtA
	lexptjA<-modelParams$precalcMatrices[[3]]$lexptjA
    }
    c(sapply(1:ncol(mSpecDist),function(s){
	.calc.mean.ouch.mv(lexpmtA[[s]],modelParams$vY0,modelParams$mPsi,modelParams$mPsi0,
	{if(is.null(modelParams$regimeTimes)){NULL}else{lexptjA[[s]]}},{if(is.null(modelParams$regimeTimes)){NULL}else{modelParams$regimes[[s]]}})
    })) 
}

.ouch.phyl.cov<-function(mTreeDist,vSpeciesDist,modelParams,vSpeciesPairs,ultrametric=FALSE){
    n<-nrow(mTreeDist)
    kY<-ncol(modelParams$A)
    mPhylCov<-matrix(NA,nrow=n*kY,ncol=n*kY)
    lCovMats<-sapply(vSpeciesPairs,function(spPair){
	s1<-(spPair-1)%/%n+1
	s2<-(spPair-1)%%n+1	
	mexpmt1A<-.calc.exptA((-1)*mTreeDist[s1,s2],modelParams$precalcMatrices[[1]])
	if(ultrametric){mexpmt2AT<-t(mexpmt1A)} ## if the tree is ultrametric then this will be the same
	else{mexpmt2AT<-t(.calc.exptA((-1)*mTreeDist[s2,s1],modelParams$precalcMatrices[[1]]))} 	
	mG<-.calc.cov.ouch.mv(vSpeciesDist[s1]-mTreeDist[s1,s2],modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])
	mexpmt1A%*%mG%*%mexpmt2AT	
    },simplify=FALSE)
    for (i in 1:length(vSpeciesPairs)){
	spPair<-vSpeciesPairs[i]
	s1<-(spPair-1)%/%n+1 
	s2<-(spPair-1)%%n+1	
	mPhylCov[((s1-1)*kY+1):(s1*kY),((s2-1)*kY+1):(s2*kY)]<-lCovMats[[i]]
	if(s1!=s2){mPhylCov[((s2-1)*kY+1):(s2*kY),((s1-1)*kY+1):(s1*kY)]<-t(lCovMats[[i]])}	
    }
    mPhylCov<-(mPhylCov+t(mPhylCov))/2
    mPhylCov
}

.mvslouch.phyl.mean<-function(mSpecDist,modelParams){
    if (is.null(modelParams$precalcMatrices[[3]])){
	lexpmtA<-sapply(mSpecDist[nrow(mSpecDist),],function(t){.calc.exptA(-t,modelParams$precalcMatrices[[1]])},simplify=FALSE) ## last row in mSpecDist will be the times of current species
        if (!(is.null(modelParams$regimeTimes))){lexptjA<-
    		    sapply(1:length(modelParams$regimeTimes),function(i,regimeTimes,vSpecDist){tjs<-regimeTimes[[i]];specT<-vSpecDist[i];sapply(tjs,function(t,specT){.calc.exptA(t-specT,modelParams$precalcMatrices[[1]])},specT=specT,simplify=FALSE)},regimeTimes=modelParams$regimeTimes,vSpecDist=mSpecDist[nrow(mSpecDist),],simplify=FALSE)}
        else{ lexptjA<-sapply(1:ncol(mSpecDist),function(i,k){list(lexpmtA[[i]],diag(1,nrow=k,ncol=k))},k=length(modelParams$vY0),simplify=FALSE)}        
    }else{
    	lexpmtA<-modelParams$precalcMatrices[[3]]$lexpmtA
	lexptjA<-modelParams$precalcMatrices[[3]]$lexptjA
    }    
    c(sapply(1:ncol(mSpecDist),function(s){
	.calc.mean.slouch.mv(lexpmtA[[s]],modelParams$precalcMatrices[[1]]$A1B,modelParams$vY0,modelParams$vX0,modelParams$mPsi,modelParams$mPsi0,
	{if(is.null(modelParams$regimeTimes)){NULL}else{lexptjA[[s]]}},{if(is.null(modelParams$regimeTimes)){NULL}else{modelParams$regimes[[s]]}})
    })) 
}


.mvslouch.phyl.cov<-function(mTreeDist,vSpeciesDist,modelParams,vSpeciesPairs,ultrametric=FALSE){
    ## the mTreeDist matrix is not symmetrical so we have to have a convention for it (tree is not ultrametric)
    ## perhaps do some simplification in case of ultrametricity
    ## if not then structure M[i,j] the time from species i to the last common ancestor of species j
    ## vSpeciesPairs is just the set of indexes but as it will always be used it is precalculated
    n<-nrow(mTreeDist)
    kY<-ncol(modelParams$A)
    kX<-ncol(modelParams$B)
    lCovMats<-sapply(vSpeciesPairs,function(spPair){
	s1<-(spPair-1)%/%n+1
	s2<-(spPair-1)%%n+1	
	mCovta<-.calc.cov.slouch.mv(vSpeciesDist[s1]-mTreeDist[s1,s2],modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])
	mexpmt1A<-.calc.exptA((-1)*mTreeDist[s1,s2],modelParams$precalcMatrices[[1]]) 
	if(ultrametric){mexpmt2AT<-t(mexpmt1A)} ## if the tree is ultrametric then this will be the same
	else{mexpmt2AT<-t(.calc.exptA((-1)*mTreeDist[s2,s1],modelParams$precalcMatrices[[1]]))} 

	mCovY1Y2<-mexpmt1A%*%mCovta[1:kY,1:kY]%*%mexpmt2AT + mexpmt1A%*%mCovta[1:kY,(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY)) + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),1:kY]%*%mexpmt2AT + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY))
	mCovY1T2<-mexpmt1A%*%mCovta[1:kY,(kY+1):(kY+kX)] + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]
	mCovT1Y2<-mCovta[(kY+1):(kY+kX),1:kY]%*%mexpmt2AT + mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY))
	mCovT1T2<-mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]
	rbind(cbind(mCovY1Y2,mCovY1T2),cbind(mCovT1Y2,mCovT1T2))
    },simplify =FALSE)
    ## Now we have a list with all of the matrix components needed -> we need to make a matrix out of it ....
    mPhylCov<-matrix(0,nrow=n*(kY+kX),ncol=n*(kY+kX))    
    for (i in 1:length(vSpeciesPairs)){
	spPair<-vSpeciesPairs[i]
	s1<-(spPair-1)%/%n+1 
	s2<-(spPair-1)%%n+1	
	mPhylCov[((s1-1)*(kY+kX)+1):(s1*(kY+kX)),((s2-1)*(kY+kX)+1):(s2*(kY+kX))]<-lCovMats[[i]]
	if(s1!=s2){mPhylCov[((s2-1)*(kY+kX)+1):(s2*(kY+kX)),((s1-1)*(kY+kX)+1):(s1*(kY+kX))]<-t(lCovMats[[i]])}
    }
    mPhylCov<-(mPhylCov+t(mPhylCov))/2
    mPhylCov
}
