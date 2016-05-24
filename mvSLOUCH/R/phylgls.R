.slouch.mv.residual.covar<-function(modelParams,B,mTreeDist,vSpeciesTimes,vSpeciesPairs,mCovPhyl=NA,invSXX=NA,mX,BFullXNA=FALSE){
    n<-length(vSpeciesTimes)
    kX<-ncol(modelParams$Sxx)
    kY<-nrow(modelParams$A)
    modelParamsB<-modelParams
    modelParamsB$B<-B
    vNAX<-which(is.na(mX))
    if ((length(vNAX)>0)&&(!BFullXNA)){
	vNAcols<-sapply(1:n,function(i,mX){if (length(which(is.na(mX[,i])))>0){i}else{-1}},mX=mX,simplify=TRUE)
	noNAcols<-which(vNAcols==-1);if (length(noNAcols)>0){vNAcols<-vNAcols[-noNAcols]}
	vNAX<-c(sapply(vNAcols,function(i,kX){((i-1)*kX+1):(i*kX)},kX=kX,simplify=TRUE))
    }
    if (is.na(mCovPhyl)){mCovPhyl<-.calc.phyl.cov(mTreeDist,vSpeciesTimes,NULL,vSpeciesPairs,"mvslouch",modelParams)}
    if (is.na(invSXX[1])){
	SXX<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX))),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX)))]
	if (length(vNAX)>0){## correct if there are missing predictor values	
	    SXX<-SXX[-vNAX,-vNAX]
	}
	invSXX<-solve(SXX)
    }   
    SYY<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY)),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY))]
    SYX<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY)),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX)))]
    SXY<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX))),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY))]

    if (length(vNAX)>0){if(nrow(SYX)>1){SYX<-SYX[,-vNAX]}else{SYX<-matrix(SYX[,-vNAX],nrow=1)};if(ncol(SXY)>1){SXY<-SXY[-vNAX,]}else{SXY<-matrix(SXY[-vNAX,],ncol=1)}}## correct if there are missing predictor values
    V<-SYY-SYX%*%invSXX%*%SXY 
    V[which(abs(V)<1e-15)]<-0
    V<-(V+t(V))/2
    list(V,mCovPhyl,invSXX)
}

.calc.vec.dist<-function(v1,v2,method="euclid",params=list(p=2)){
    dist<-0
    if (method=="euclid"){dist<-(sum((v1-v2)^params$p))^(1/params$p)}
    dist
}

.solve.reg<-function(D,Y,V,UnknownIntercept=FALSE,method="psinv"){
## function uses Cholesky decomposition and then backsubstatution to solve the problem
    vNAY<-which(is.na(Y)) ## remove missing value rows
    if (length(vNAY)>0){Y<-Y[-vNAY];if(ncol(D)>1){D<-D[-vNAY,]}else{D<-matrix(D[-vNAY,],ncol=1)};V<-V[-vNAY,-vNAY]}
    sol<-matrix(NA,nrow=ncol(D),ncol=1)    
    if (ncol(D)>0){
	if (method == "chol"){
	    tryCatch({
		    invV<-solve(V)
		    W<-t(D)%*%invV%*%D
	    	    R<-chol(W)
		    z<-solve(t(R))%*%t(D)%*%Y
		    sol<-solve(R)%*%z
		},
		error=function(e){print(paste("Cholesky ",e))}
	    )
	}    
	if (method == "psinv"){
	    tryCatch({
		    vP<-chol(V) ## make the transormation to have independent residuals, as described in
		    Pt<-solve(t(vP)) ## Box Jenkins, Time Series Analysis Holden-Day San Francisco 1970 p265-267
	    	    pD<-Pt%*%D
		    pY<-Pt%*%Y		
		    psinvpD<-pseudoinverse(pD) ## pseudoinverse comes from corpcor library
		    sol<-psinvpD%*%pY 
		},
		error=function(e){print(paste("pseudo inverse",e))}
	    )
	}
    }
    sol
}

.slouch.mv.iGLS<-function(lPrecalculates,Y,modelParams,designToEstim,mX,UnknownIntercept=FALSE,X0=NULL,tol=0.0001,maxIter=50,ShouldPrint=TRUE){
    lDzeta<-modelParams$precalcMatrices[[4]]$lDzeta
    lDzetaKappa<-modelParams$precalcMatrices[[4]]$lDzetaKappa
    lexpmtA<-modelParams$precalcMatrices[[3]]$lexpmtA
    lexptjA<-modelParams$precalcMatrices[[3]]$lexptjA
    intercept<-modelParams$precalcMatrices[[5]]$intercept
    invmAncestorTimes<-lPrecalculates$invmAncestorTimes
    mSpecDist<-lPrecalculates$mSpecDist
    mTreeDist<-lPrecalculates$mTreeDist
    mAncestorTimes<-lPrecalculates$mAncestorTimes
    vSpeciesPairs<-lPrecalculates$vSpeciesPairs
    
    kY<-nrow(modelParams$A)
    kX<-ncol(modelParams$Syx)
    
    vNAcols<-c()
    if (!designToEstim$FullNAYX){
	vNAcols<-sapply(1:ncol(mX),function(i,mX){if (length(which(is.na(mX[,i])))>0){i}else{-1}},mX=mX,simplify=TRUE)
	noNAXcols<-which(vNAcols==-1);if (length(noNAXcols)>0){vNAcols<-vNAcols[-noNAXcols]}	
    }
    if (!designToEstim$FullNAY){
	mY<-matrix(Y,ncol=ncol(mX))
	vYNAcols<-sapply(1:ncol(mX),function(i,mY){if (length(which(is.na(mY[,i])))>0){i}else{-1}},mY=mY,simplify=TRUE)
	noNAYcols<-which(vYNAcols==-1);if (length(noNAYcols)>0){vYNAcols<-vYNAcols[-noNAYcols]}	
	vNAcols<-union(vYNAcols,vNAcols)
    }
    if (length(vNAcols)>0){Y[c(sapply(vNAcols,function(i,kY){((i-1)*kY+1):(i*kY)},kY=kY,simplify=TRUE))]<-NA}
    
    lDesign<-.slouch.mv.design.matrix.YcT(mSpecDist,mTreeDist,mAncestorTimes,modelParams,designToEstim,mX,lDzeta,lDzetaKappa,lexpmtA,lexptjA,X0,invmAncestorTimes)
    Bcols<-lDesign$Bcol:(lDesign$Bcol+kX*kY-1)    
    modelParams$intercept<-intercept
    Y<-Y-intercept
    if (designToEstim$B){	
	if (UnknownIntercept){lDesign$colB<-lDesign$colB+kY}	
	mCovPhyl<-modelParams$mCovPhyl
	invSXX<-modelParams$invSXX
	V<-diag(1,nrow=nrow(lDesign$D),ncol=nrow(lDesign$D)) ## we assume first that everything is independent
	estimParams<-.solve.reg(D=lDesign$D,Y=Y,V=V,UnknownIntercept)	
	estimParams1<-estimParams+tol*100000000000
	iter<-1    	    
	## ----------------------------------------------------------------------
	if (designToEstim$iRegLin){	
	    while((length(which(is.na(estimParams)))==0)&&(.calc.vec.dist(estimParams,estimParams1)>tol[1])&&(iter<=maxIter)){
		if (ShouldPrint){print(paste("Running iteration ",iter,"of iterated GLS"))}
		estimParams1<-estimParams
		B<-matrix(estimParams[Bcols,1],nrow=kY,ncol=kX,byrow=TRUE) 
		modelParams$B<-B
		modelParams$precalcMatrices[[1]]$A1B<-modelParams$precalcMatrices[[1]]$invA%*%B 
		lCovars<-.slouch.mv.residual.covar(modelParams,B,mTreeDist,mSpecDist[nrow(mSpecDist),],vSpeciesPairs,NA,invSXX,mX,designToEstim$BFullXNA)	
		V<-lCovars[[1]]
		mCovPhyl<-lCovars[[2]]
		invSXX<-lCovars[[3]]
		estimParams<-.solve.reg(D=lDesign$D,Y=Y,V=V,UnknownIntercept)	
		if (is.na(estimParams[1])){estimParams<-estimParams1} ## we had an error -> get previous estimate and get out of this point
		iter<-iter+1	    	    
    	    }
        }
        else{	    
	    par0<-estimParams[,1]
	    if (designToEstim$optimMethod=="optim"){	    
	    	tryCatch({
	    	    OptimRes<-optim(par=par0,fn=.min.RSS.BmPsi,D=lDesign$D,Y=Y,estimParams=estimParams,Bcols=Bcols,modelParams=modelParams,mTreeDist=mTreeDist,mSpecDist=mSpecDist,vSpeciesPairs=vSpeciesPairs,invSXX=invSXX,kY=kY,kX=kX)
		},error=function(e){print(e)})
		estimParams<-matrix(OptimRes$par,ncol=1)    
	    }
	    if (designToEstim$optimMethod=="nlminb"){
	    	tryCatch({
	    	    OptimRes<-nlminb(start=par0,objective=.min.RSS.BmPsi,D=lDesign$D,Y=Y,estimParams=estimParams,Bcols=Bcols,modelParams=modelParams,mTreeDist=mTreeDist,mSpecDist=mSpecDist,vSpeciesPairs=vSpeciesPairs,invSXX=invSXX,kY=kY,kX=kX)
		},error=function(e){print(e)})
		estimParams<-matrix(OptimRes$par,ncol=1)
	    }
        }
        B<-matrix(estimParams[lDesign$Bcol:(lDesign$Bcol+kX*kY-1),1],nrow=kY,ncol=kX,byrow=TRUE)
	modelParams$B<-B
	modelParams$precalcMatrices[[1]]$A1B<-modelParams$precalcMatrices[[1]]$invA%*%B
	lCovars<-.slouch.mv.residual.covar(modelParams,B,mTreeDist,mSpecDist[nrow(mSpecDist),],vSpeciesPairs,NA,invSXX,mX,designToEstim$BFullXNA)	
	mCovPhyl<-lCovars[[2]]
	invSXX<-lCovars[[3]]	
	if (ShouldPrint){
	    if ((length(which(is.na(estimParams)))!=0)||(.calc.vec.dist(estimParams,estimParams1)>=tol[1])){print("Iterated GLS did not reach convergence at point")}
	    else{print("Iterated GLS converged")}
	}
    }
    else{
        n<-ncol(mSpecDist)
        if (is.na(modelParams$mCovPhyl[1])){mCovPhyl<-.calc.phyl.cov(mTreeDist,mSpecDist[nrow(mSpecDist),],NULL,vSpeciesPairs,"mvslouch",modelParams)}
	else{mCovPhyl<-modelParams$mCovPhyl}
	SYY<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY)),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY))]
	SYX<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY)),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX)))]
        SXY<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX))),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY))]
	if (is.na(modelParams$invSXX[1])){
	    SXX<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX))),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX)))]
	    invSXX<-solve(SXX)
	}
	else{invSXX<-modelParams$invSXX}
	V<-SYY-SYX%*%invSXX%*%SXY 
	V[which(abs(V)<1e-15)]<-0
        estimParams<-.solve.reg(D=lDesign$D,Y=Y,V=V,UnknownIntercept)
    }
    estimParams[which(abs(estimParams)<1e-15)]<-0
    if (UnknownIntercept){modelParams$EstimedIntercept<-estimParams[1:kY,1];CurrPos<-kY+1}else{CurrPos<-1}       
    if (designToEstim$y0 && !designToEstim$y0AncState){modelParams$vY0<-estimParams[1:kY,1];CurrPos<-kY+1}
    
    ## we get y0 last as we might be mapping it back to optimal ancestral state
    ## get psi
    if (designToEstim$psi){
	modelParams$mPsi<-matrix(estimParams[CurrPos:(CurrPos+length(modelParams$regimeTypes)*kY-1),1],nrow=kY,ncol=length(modelParams$regimeTypes),byrow=FALSE)
	CurrPos<-CurrPos+length(modelParams$regimeTypes)*kY
    }
    ## get psi0
    if (designToEstim$psi0){
	modelParams$mPsi0<-matrix(estimParams[CurrPos:(CurrPos+kY-1),1],nrow=kY,ncol=1)
	CurrPos<-CurrPos+kY
    }
    ## get B
    if (designToEstim$B){
	modelParams$B<-matrix(estimParams[Bcols,1],nrow=kY,ncol=kX,byrow=TRUE) 
	modelParams$precalcMatrices[[1]]$A1B<-modelParams$precalcMatrices[[1]]$invA%*%modelParams$B
	CurrPos<-CurrPos+kX*kY
    }
    
    ## get X0    
    if (designToEstim$X0){modelParams$vX0<-estimParams[CurrPos:(CurrPos+kX-1),1]}

    ## get y0
    if (designToEstim$y0 && designToEstim$y0AncState){
    	modelParams$vY0<-modelParams$mPsi[,designToEstim$y0Regime]+modelParams$mPsi0
        if (!designToEstim$y0OnlyFixed){modelParams$vY0<-modelParams$vY0-modelParams$precalcMatrices[[1]]$A1B%*%modelParams$vX0}
    }    
    modelParams$mCovPhyl<-mCovPhyl
    modelParams$invSXX<-invSXX
    vNAY<-which(is.na(Y));if (length(vNAY)>0){if(ncol(lDesign$D)>1){lDesign$D<-lDesign$D[-vNAY,]}else{lDesign$D<-matrix(lDesign$D[-vNAY,],ncol=1)};V<-V[-vNAY,-vNAY]}

    modelParams$regressCovar<-matrix(0,nrow=ncol(lDesign$D),ncol=ncol(lDesign$D))
    if(ncol(lDesign$D)>0){modelParams$regressCovar<-pseudoinverse(t(lDesign$D)%*%pseudoinverse(V)%*%lDesign$D)}
    modelParams    
}


.slouch.mv.design.matrix.YcT<-function(mSpecDist,mTreeDist,mAncestorTimes,modelParams,designToEstim,mX,lDzeta=NULL,lDzetaKappa=NULL,lexpmtA=NULL,lexptjA=NULL,X0=NULL,invmAncestorTimes=NULL){
    invA<-modelParams$precalcMatrices[[1]]$invA
    A<-modelParams$A
    mCovXX<-modelParams$precalcMatrices[[2]]$S22 
    mCovTX<-modelParams$precalcMatrices[[2]]$S12 
    kY<-nrow(A) ## number of traits
    kX<-ncol(mCovTX)
    n<-nrow(mTreeDist) ## number of species
    iNumToEstim<-0
    if (designToEstim$y0 && !designToEstim$y0AncState){iNumToEstim<-iNumToEstim+kY} ## original trait
    if (designToEstim$psi){ iNumToEstim<-iNumToEstim+length(modelParams$regimeTypes)*kY }
    if (designToEstim$psi0){ iNumToEstim<-iNumToEstim+kY }
    if (designToEstim$B){ iNumToEstim<-iNumToEstim+kY*kX }
    if (designToEstim$X0){iNumToEstim<-iNumToEstim+kX} ## originial predictor
    if (designToEstim$BX0){iNumToEstim<-iNumToEstim+kY} ## vector BX0 
    D<-matrix(0,nrow=n*kY,ncol=iNumToEstim) 
    currXcol<-1
	
    ## setup Y0 column ---------------------------------------------------------------------------
    if (designToEstim$y0 && !designToEstim$y0AncState){ 
	D[,1:kY]<-sapply(lexpmtA,function(x){x},simplify=TRUE)
	currXcol<-currXcol+kY
    }
    #---------------------------------------------------------------------------------------------
    ## seutp psi ---------------------------------------------------------------------------------
    ## regime names assumed to be numbers 1 .. max and according to this constructed
    if (designToEstim$psi){ 
	for (i in 1:n){## for each species
	    mRegTmp<-matrix(0,ncol=length(modelParams$regimeTypes)*kY,nrow=kY)
	    for (j in 1:length(modelParams$regimes[[i]])){
	    	jType<-modelParams$regimes[[i]][j]
		mRegTmp[,((jType-1)*kY+1):(jType*kY)]<-mRegTmp[,((jType-1)*kY+1):(jType*kY)]+(lexptjA[[i]][[j+1]]-lexptjA[[i]][[j]])
	    }
	    if (designToEstim$y0 && designToEstim$y0AncState){ 
    		jType=designToEstim$y0Regime
    		mRegTmp[,((jType-1)*kY+1):(jType*kY)]<-mRegTmp[,((jType-1)*kY+1):(jType*kY)]+lexpmtA[[i]]
    	    }
	    D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+length(modelParams$regimeTypes)*kY-1)]<-mRegTmp
	}
	currXcol<-currXcol+length(modelParams$regimeTypes)*kY    
    }
    #---------------------------------------------------------------------------------------------    
    ## setup Psi0 column ---------------------------------------------------------------------------
    if (designToEstim$psi0){ 	
	if (designToEstim$y0 && designToEstim$y0AncState){lPsi0Coeffs<-sapply(1:n,function(i,kY){diag(1,kY,kY)},simplify=FALSE)}
	else{lPsi0Coeffs<-sapply(lexpmtA,function(x){diag(1,ncol(x),nrow(x))-x},simplify=FALSE)}
    	for (i in 1:n){## for each species
    	    D[((i-1)*kY+1):(i*kY),currXcol:(currXcol+kY-1)]<-lPsi0Coeffs[[i]]
    	}
	currXcol<-currXcol+kY
    }
    #---------------------------------------------------------------------------------------------

    ## seutp B -----------------------------------------------------------------------------------
    if (designToEstim$B){
	for (i in 1:n){## for each species
	    if (designToEstim$UseX0){D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]<-lDzeta[[i]]%x%matrix(X0,nrow=1)}
	    if (designToEstim$SimpReg){
		vToAdd<-rep(0,kX*kY)
		viXmX0<-mX[,i]-X0
		vNAX<-which(is.na(viXmX0))
		if (length(vNAX)>0){
		    if ((length(vNAX)<kX)&&designToEstim$BFullXNA){
			viXmX0<-viXmX0[-vNAX]
			tmpimDzetaKappa<-(lDzetaKappa[[i]]%*%solve(modelParams$precalcMatrices[[2]]$S22))
			if(nrow(tmpimDzetaKappa)>1){imDzetaKappa<-tmpimDzetaKappa[,-vNAX]}else{imDzetaKappa<-matrix(tmpimDzetaKappa[,-vNAX],nrow=1)}
			vToAdd<-imDzetaKappa%x%matrix(viXmX0,nrow=1) 
		    }
		}else{vToAdd<-lDzetaKappa[[i]]%x%matrix(viXmX0,nrow=1)}		
		D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]<-D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]+vToAdd
	    }
	    else{D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]<-D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]+lDzetaKappa[[i]]}
	    if (designToEstim$y0 && designToEstim$y0AncState && !designToEstim$y0OnlyFixed && designToEstim$UseX0){D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]<-D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]-(lexpmtA[[i]]%*%invA)%x%matrix(X0,nrow=1)}    
	}
	Bcol<-currXcol
	currXcol<-currXcol+kY*kX    
    }else{Bcol<-NA}
    #---------------------------------------------------------------------------------------------    		
    ## setup X0 ----------------------------------------------------------------------------------
    if (designToEstim$X0){
	for (i in 1:n){## for each species
	    invCovXX<-modelParams$precalcMatrices[[2]]$invS22
	    X0coeff<-sapply(1:n,function(j,i,precalcsA,kY,mTreeDist,mAncestorTimes,invmAncestorTimes){
		.calc.exptA(-mTreeDist[i,j],precalcsA)%*%(.calc.exptA(mAncestorTimes[i,j],precalcsA)-diag(1,kY,kY))*sum(invmAncestorTimes[j,])
	    },i=i,precalcsA=modelParams$precalcMatrices[[1]],kY=kY,mTreeDist=mTreeDist,mAncestorTimes=mAncestorTimes,invmAncestorTimes=invmAncestorTimes,simplify=FALSE)
	    X0coeff<-Reduce('+',X0coeff,matrix(0,nrow=kY,ncol=kX))%*%mCovTX%*%invCovXX
	    D[((i-1)*kY+1):(i*kY),currXcol:(currXcol+kX-1)]<-(-1)*X0coeff	    
	}
	X0colStart<-currXcol
	currXcol<-currXcol+kX
    }
    ## -------------------------------------------------------------------------------------------
    ## setup B*X0 --------------------------------------------------------------------------------    
    if (designToEstim$BX0){
	for (i in 1:n){## for each species
	    D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY-1)]<-lDzeta[[i]]
	    D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY-1)]<-D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY-1)]-
	    Reduce('+',sapply(1:n,function(j,i,precalcsA,invA,mTreeDist,mAncestorTimes,invmAncestorTimes,kY){
		    (.calc.exptA(-mTreeDist[i,j],precalcsA)%*%invA%*%(.calc.exptA(mAncestorTimes[i,j],precalcsA)-diag(1,kY,kY))-diag(mAncestorTimes[i,j],kY,kY)
		    )%*%invA*sum(invmAncestorTimes[j,])}
		,i=i,precalcsA=modelParams$precalcMatrices[[1]],invA=invA,mTreeDist=mTreeDist,mAncestorTimes=mAncestorTimes,invmAncestorTimes=invmAncestorTimes,kY=kY
	    ,simplify=FALSE),rep(0,kY))	
	}
    }    
    ## -------------------------------------------------------------------------------------------
    D[which(abs(D)<1e-15)]<-0
    list("D"=D,"Bcol"=Bcol)
}

.calc.intercept<-function(){

}

.interpret.intercept<-function(){
}

.ouch.GLS<-function(lPrecalculates,Y,modelParams,designToEstim,tol=0.0001,UnknownIntercept=FALSE,intercept=NA){
    lexpmtA<-modelParams$precalcMatrices[[3]]$lexpmtA
    lexptjA<-modelParams$precalcMatrices[[3]]$lexptjA
    mSpecDist<-lPrecalculates$mSpecDist
    mTreeDist<-lPrecalculates$mTreeDist
    mAncestorTimes<-lPrecalculates$mAncestorTimes
    vSpeciesPairs<-lPrecalculates$vSpeciesPairs

    kY<-nrow(modelParams$A)

    if (!designToEstim$FullNAY){
	mY<-matrix(Y,ncol=ncol(mTreeDist))
	vYNAcols<-sapply(1:ncol(mTreeDist),function(i,mY){if (length(which(is.na(mY[,i])))>0){i}else{-1}},mY=mY,simplify=TRUE)
	noNAYcols<-which(vYNAcols==-1);if (length(noNAYcols)>0){vYNAcols<-vYNAcols[-noNAYcols]}	
    	Y[c(sapply(vYNAcols,function(i,kY){((i-1)*kY+1):(i*kY)},kY=kY,simplify=TRUE))]<-NA
    }

    if (!is.na(intercept)){Y<-Y-intercept}

    lDesign<-.ouch.design.matrix(mSpecDist,mTreeDist,modelParams,designToEstim,lexpmtA,lexptjA)
    if (!is.na(modelParams$mCovPhyl)){V<-modelParams$mCovPhyl}else{V<-.calc.phyl.cov(mTreeDist,mSpecDist[nrow(mSpecDist),],NULL,vSpeciesPairs,"ouch",modelParams)}
    estimParams<-.solve.reg(D=lDesign$D,Y=Y,V=V,UnknownIntercept)

    estimParams[which(abs(estimParams)<1e-15)]<-0
    if (UnknownIntercept){modelParams$EstimedIntercept<-estimParams[1:kY,1];CurrPos<-kY+1}else{CurrPos<-1}   
    if (designToEstim$y0 && !designToEstim$y0AncState){modelParams$vY0<-estimParams[1:kY,1];CurrPos<-kY+1}    

    if (designToEstim$psi){
        modelParams$mPsi<-matrix(estimParams[CurrPos:(CurrPos+length(modelParams$regimeTypes)*kY-1),1],nrow=kY,ncol=length(modelParams$regimeTypes),byrow=FALSE)
        CurrPos<-CurrPos+length(modelParams$regimeTypes)*kY
    }
    if (designToEstim$psi0){
        modelParams$mPsi0<-matrix(estimParams[CurrPos:(CurrPos+kY-1),1],nrow=kY,ncol=1)
        CurrPos<-CurrPos+kY
    }
    if (designToEstim$y0 && designToEstim$y0AncState){modelParams$vY0<-matrix(modelParams$mPsi[,designToEstim$y0Regime]+modelParams$mPsi0,ncol=1,nrow=kY)}    
    modelParams$mCovPhyl<-V    
    vNAY<-which(is.na(Y));if (length(vNAY)>0){if(ncol(lDesign$D)>1){lDesign$D<-lDesign$D[-vNAY,]}else{lDesign$D<-matrix(lDesign$D[-vNAY,],ncol=1)};V<-V[-vNAY,-vNAY]}
    modelParams$regressCovar<-matrix(0,nrow=ncol(lDesign$D),ncol=ncol(lDesign$D))
    if(ncol(lDesign$D)>0){modelParams$regressCovar<-pseudoinverse(t(lDesign$D)%*%pseudoinverse(V)%*%lDesign$D)}
    
    modelParams
}

.ouch.design.matrix<-function(mSpecDist,mTreeDist,modelParams,designToEstim,lexpmtA=NULL,lexptjA=NULL){
## order of columns in mX, has to be the same as order of species in phylogeny
    invA<-modelParams$precalcMatrices[[1]]$invA
    A<-modelParams$A
    kY<-nrow(A) ## number of traits
    n<-nrow(mTreeDist) ## number of species

    iNumToEstim<-0
    if (designToEstim$y0 && !designToEstim$y0AncState){iNumToEstim<-iNumToEstim+kY} ## original trait
    if (designToEstim$psi){ iNumToEstim<-iNumToEstim+length(modelParams$regimeTypes)*kY }
    if (designToEstim$psi0){ iNumToEstim<-iNumToEstim+kY }
    D<-matrix(0,nrow=n*kY,ncol=iNumToEstim)

    currXcol<-1
    
    ## setup Y0 column ---------------------------------------------------------------------------
    if (designToEstim$y0 && !designToEstim$y0AncState){ 
	D[,1:kY]<-sapply(lexpmtA,function(x){x},simplify=TRUE) 
	currXcol<-currXcol+kY
    }
    #---------------------------------------------------------------------------------------------
    ## setup psi ---------------------------------------------------------------------------------
    ## regime names assumed to be numbers 1 .. max and according to this constructed
    if (designToEstim$psi){ 
	for (i in 1:n){## for each species
	    mRegTmp<-matrix(0,ncol=length(modelParams$regimeTypes)*kY,nrow=kY)
	    for (j in 1:length(modelParams$regimes[[i]])){
		jType<-modelParams$regimes[[i]][j]
		mRegTmp[,((jType-1)*kY+1):(jType*kY)]<-mRegTmp[,((jType-1)*kY+1):(jType*kY)]+(lexptjA[[i]][[j+1]]-lexptjA[[i]][[j]])
	    }
	    if (designToEstim$y0 && designToEstim$y0AncState){ 
    		jType=designToEstim$y0Regime
    		mRegTmp[,((jType-1)*kY+1):(jType*kY)]<-mRegTmp[,((jType-1)*kY+1):(jType*kY)]+lexpmtA[[i]]
    	    }
	    D[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+length(modelParams$regimeTypes)*kY-1)]<-mRegTmp
	}
	currXcol<-currXcol+length(modelParams$regimeTypes)*kY
    }
    #---------------------------------------------------------------------------------------------    
    ## setup Psi0 column ---------------------------------------------------------------------------
    if (designToEstim$psi0){ 
	D[,currXcol:(currXcol+kY-1)]<-sapply(lexpmtA,function(x){diag(1,ncol(x),nrow(x))-x},simplify=TRUE)
	currXcol<-currXcol+kY
    }
    #---------------------------------------------------------------------------------------------

    D[which(abs(D)<1e-15)]<-0
    list("D"=D)
}

