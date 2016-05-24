.maxlik.estim<-function(dfData,dfData.Merror,PhylTree,EvolModel,EstimationParams,regimeTimes=NULL,regimes=NULL,regimeTypes=NULL,method="igls",tol=0.0001,maxIter=50,bShouldPrint=FALSE,maxTries=10,minLogLik=-Inf){
    MaxLikEstim<-NA
    if (is.null(EstimationParams$vVars)){EstimationParams$vVars<-NULL}
    if (is.null(EstimationParams$conditional)){EstimationParams$conditional<-FALSE}
    if (is.null(EstimationParams$Atype)){EstimationParams$Atype<-"DecomposableReal"}
    if (is.null(EstimationParams$Btype)){EstimationParams$Btype<-"Any"}
    if (is.null(EstimationParams$mPsitype)){EstimationParams$mPsitype<-"Global"}
    if (is.null(EstimationParams$Syytype)){EstimationParams$Syytype<-"Symmetric"}
    if (is.null(EstimationParams$Sxxtype)){EstimationParams$Sxxtype<-"Symmetric"}
    
    lPrecalculates<-.calculate.Tree.dists(PhylTree,UserTermLabels=EstimationParams$TerminalLabels)
    if (is.null(regimeTimes)){
	vSpT<-lPrecalculates$mSpecDist[1,]
	names(vSpT)<-c()
	regimeTimes<-sapply(vSpT,function(x){c(0,x)},simplify=FALSE)
    }
    if ((is.null(regimes))||(is.null(regimeTypes))){
	regimes<-as.list(rep("1",PhylTree@nterm))
	regimeTypes<-c("1")
    }

    modelParams<-vector("list",4)
    names(modelParams)<-c("regimeTimes","regimes","regimeTypes","Merror")
    modelParams$regimeTimes<-regimeTimes
    modelParams$regimes<-regimes
    modelParams$regimeTypes<-regimeTypes		
    
    ## dfData is assumed to be in the following format 
    ## number of rows is the number of species, first columns 1:kY response, columns (kY+1):(kY+kX) predictor (for mvslouch)
    ## needs to be done with t as c(M[,cols]) vectorizes M[,cols] by column
    data<-c(t(dfData))
    modelParams$Merror<-matrix(0,length(data),length(data))
    ## we assume at the moment that the measurement errors are independent of each other
    ##if (!is.null(dfData.Merror)){modelParams$Merror<-diag(c(t(dfData.Merror)),length(data),length(data))}
    if (!is.null(dfData.Merror)){modelParams$Merror<-.createCovariancematrix(dfData.Merror,nrow(dfData),ncol(dfData),NULL,"measurement error")}
    ##print(modelParams$Merror)
    ##print(dim(modelParams$Merror))
    vY<-c(t(dfData)[1:EstimationParams$kY,])
    if ((method=="maxlik")&&(EvolModel=="bm")){
	dfData<-as.matrix(dfData)
	dfData<-as.data.frame(rbind(matrix(NA,ncol=ncol(dfData),nrow=PhylTree@nnodes-PhylTree@nterm),dfData))
	rownames(dfData)<-1:nrow(dfData)                
	MaxLikEstim<-vector("list",3)
	names(MaxLikEstim)<-c("BrownResult","ParamsInModel","ParamSummary")
	MaxLikEstim$BrownResult<-.bm.estim(dfData,PhylTree,modelParams$Merror)
	MaxLikEstim$ParamsInModel<-list("Sxx"=MaxLikEstim$BrownResult$Sxx,"vX0"=MaxLikEstim$BrownResult$vX0)
	if (EstimationParams$calcCI){
             if(bShouldPrint){print("Calculating confidence intervals can take very long time")}
	     MaxLikEstim$ParamSummary<-.Params.summary(MaxLikEstim$ParamsInModel,"bm",NULL,dfData[(PhylTree@nnodes-PhylTree@nterm+1):PhylTree@nnodes,],NULL,MaxLikEstim$BrownResult$LogLik,PhylTree@nterm,NULL,MaxLikEstim$BrownResult$RSS,lPrecalculates=lPrecalculates,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik)
        }else{MaxLikEstim$ParamSummary<-.Params.summary(MaxLikEstim$ParamsInModel,"bm",NULL,dfData[(PhylTree@nnodes-PhylTree@nterm+1):PhylTree@nnodes,],NULL,MaxLikEstim$BrownResult$LogLik,PhylTree@nterm,NULL,MaxLikEstim$BrownResult$RSS,lPrecalculates=list(tree.height=lPrecalculates$tree.height))}
    }    
    if ((method=="grid")||(method=="gridigls")||(method=="gridgls")){## setup grid search
	if (is.null(EstimationParams$mGrid)){mGrid<-.GenerateModelGrid(EstimationParams)}
	else{mGrid<-EstimationParams$mGrid}
	mGrid<-cbind(1:nrow(mGrid),mGrid) ## for diagnostic purposes ....
	if (bShouldPrint){print(paste("Size of grid : ",nrow(mGrid),sep=""))}
	mLogLikSurface<-NA
	if (method=="grid"){
	    mLogLikSurface<-cbind(mGrid,"LogLik"=apply(mGrid,1,
		    function(x){		    
			LogLik<-.calcLogLikParametrized(x[-1],data,vY,PhylTree,EvolModel,lPrecalculates,EstimationParams,modelParams,tol,maxIter,bShouldPrint)
			if (x[1]%%1000==1){names(x)<-c();print(c(x,LogLik),digits=2)} ##print every 1000 
			LogLik
		    }
		)
	    )
	}
	if (((method=="gridigls")&&(EvolModel=="mvslouch"))||((method=="gridgls")&&(EvolModel=="ouch"))){
	    EstimationParams<-.beginEstimationParams(EvolModel,EstimationParams,dfData,PhylTree,modelParams$MError)		   
	    lLogLikSurfaces<-apply(mGrid,1,function(x,EstimationParams,lPrecalculates,modelParams,vY,tol,maxIter,data)
	    {
	    	modelParams<-.par.transform(x[-1],EstimationParams,modelParams)
		modelParams$precalcMatrices<-.calcAprecalcMatrices(EvolModel,modelParams,lPrecalculates,EstimationParams)
		modelParams<-.GridEvaluatePointA(EvolModel,vY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint)			    
		LogLik<-.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional)
		if (x[1]%%30==1){names(x)<-c();print(c(x,LogLik),digits=2)} ##print every 1000 
		LogLik		
	    },EstimationParams,lPrecalculates,modelParams,vY,tol,maxIter,data
	    )
	    mLogLikSurface<-cbind(mGrid,"LogLik"=lLogLikSurfaces)
	}	
	MaxLik<-vector("list",3)
	names(MaxLik)<-c("GridResult","ParamsInModel","ParamSummary")
	vFiniteLiks<-intersect(intersect(which(!(is.infinite(mLogLikSurface[,"LogLik"]))),which(!(is.nan(mLogLikSurface[,"LogLik"])))),which(!(is.na(mLogLikSurface[,"LogLik"])))) 
	MaxLik$GridResult<-matrix(mLogLikSurface[which(mLogLikSurface[,"LogLik"]==max(mLogLikSurface[vFiniteLiks,"LogLik"])),],nrow=length(which(mLogLikSurface[,"LogLik"]==max(mLogLikSurface[vFiniteLiks,"LogLik"]))))
	colnames(MaxLik$GridResult)<-colnames(mLogLikSurface)
	MaxLik$ParamsInModel<-vector("list",nrow(MaxLik$GridResult))
	MaxLik$ParamSummary<-vector("list",nrow(MaxLik$GridResult))
	for (i in 1:nrow(MaxLik$GridResult)){
	    MaxLik$ParamsInModel[[i]]<-.par.transform(MaxLik$GridResult[i,-c(1,ncol(MaxLik$GridResults)-1,ncol(MaxLik$GridResults))],EstimationParams,modelParams)
	    RSS<-NA
	    if ((method=="gridigls")||(method=="gridgls")){
		MaxLik$ParamsInModel[[i]]<-.EvaluatePoint(EvolModel,data,vY,MaxLik$ParamsInModel[[i]],lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,TRUE,NULL,FALSE)$modelParams			
		RSS<-.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=MaxLik$ParamsInModel[[i]],vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,TRUE)
	    }		    	
	     if (EstimationParams$calcCI){
                if(bShouldPrint){print("Calculating confidence intervals can take very long time")}
	        MaxLik$ParamSummary[[i]]<-.Params.summary(MaxLik$ParamsInModel[[i]],EvolModel,EstimationParams$designToEstim,data,1,MaxLik$GridResult[i,"LogLik"],nrow(dfData),length(MaxLik$GridResult[i,])-2,RSS,lPrecalculates=lPrecalculates,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik)
    	    }else{MaxLik$ParamSummary[[i]]<-.Params.summary(MaxLik$ParamsInModel[[i]],EvolModel,EstimationParams$designToEstim,data,1,MaxLik$GridResult[i,"LogLik"],nrow(dfData),length(MaxLik$GridResult[i,])-2,RSS,lPrecalculates=list(tree.height=lPrecalculates$tree.height))}
	    MaxLik$ParamsInModel[[i]]<-.cleanUpModelParams(MaxLik$ParamsInModel[[i]])
	}
	MaxLikEstim<-list("MaxLik"=MaxLik,"LogLikSurface"=mLogLikSurface)
    }
    if (method=="glsgc"){
        MaxLikEstim<-.glsgc.estim(dfData=dfData,data=data,EvolModel=EvolModel,PhylTree=PhylTree,EstimationParams=EstimationParams,modelParams=modelParams,lPrecalculates=lPrecalculates,tol=tol,maxIter=maxIter,bShouldPrint=bShouldPrint,maxTries=maxTries,minLogLik=minLogLik)
    }
    
    if ((method=="EM")&&(EvolModel=="mvslouch")){
	  print("Method not yet implemented")
	  MaxLikEstim<-NULL
##        MaxLikEstim<-MVslouchEM.estim(dfData=dfData,data=data,PhylTree=PhylTree,EstimationParams=EstimationParams,modelParams=modelParams,lPrecalculates=lPrecalculates,tol=tol,maxIter=maxIter,bShouldPrint=bShouldPrint)
    }
    
    if ((method=="MCEM")&&(EvolModel=="mvslouch")){
	  print("Method not yet implemented")
	  MaxLikEstim<-NULL
##        MaxLikEstim<-MVslouchMCEM.estim(dfData=dfData,data=data,PhylTree=PhylTree,EstimationParams=EstimationParams,modelParams=modelParams,lPrecalculates=lPrecalculates,tol=tol,maxIter=maxIter,bShouldPrint=bShouldPrint)
    }
    
    if (((method=="gls")&&(EvolModel=="ouch"))||((method=="igls")&&(EvolModel=="mvslouch"))){
	if(is.null(EstimationParams$StartPoint)){print("No starting position for optimization. Cannot continue")}
	else{
	    EstimationParams<-.beginEstimationParams(EvolModel,EstimationParams,dfData,PhylTree,modelParams$Merror)		    
	    vParFound<-rep(NA,length(EstimationParams$StartPoint))	
	    names(vParFound)<-names(EstimationParams$StartPoint)
	    MaximEstim<-NA
	    par0<-EstimationParams$StartPoint
	    if (EstimationParams$maximMethod=="optim"){## setup likelihood optim	    		
		MaximEstim<-optim(
		    par=par0,
		    fn=function(par,data,PhylTree,EvolModel,lPrecalculates,EstimationParams,modelParams,minLogLik,parNames){
			(-1)*.calcLogLikParametrized(par,data,vY,PhylTree,EvolModel,lPrecalculates,EstimationParams,modelParams,TRUE,tol,maxIter,bShouldPrint,minLogLik=minLogLik,parNames=parNames)
		    },data=data,PhylTree=PhylTree,EvolModel=EvolModel,lPrecalculates=lPrecalculates,EstimationParams=EstimationParams,modelParams=modelParams,minLogLik=minLogLik,parNames=names(par0)
		)	    	    
		vParFound<-MaximEstim$par
	    }
	    if(EstimationParams$maximMethod=="nlminb"){
	        MaximEstim<-nlminb(
            	    start=par0,
            	    objective=function(par,data,PhylTree,EvolModel,lPrecalculates,EstimationParams,modelParams,minLogLik,parNames){
			(-1)*.calcLogLikParametrized(par,data,vY,PhylTree,EvolModel,lPrecalculates,EstimationParams,modelParams,TRUE,tol,maxIter,bShouldPrint,minLogLik=minLogLik,parNames=parNames)
		    },data=data,PhylTree=PhylTree,EvolModel=EvolModel,lPrecalculates=lPrecalculates,EstimationParams=EstimationParams,modelParams=modelParams,minLogLik=minLogLik,parNames=names(par0)
        	)
        	vParFound<-MaximEstim$par
	    }
	    names(vParFound)<-names(EstimationParams$StartPoint)
	    tmpModelParams<-.EvaluatePoint(EvolModel,data,vY,.par.transform(vParFound,EstimationParams,modelParams),lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
	    RSS<-.calc.phyl.LogLik.traits(data,lPrecalculates=lPrecalculates,EvolModel,modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,TRUE)
	    if (EstimationParams$calcCI){
        	if(bShouldPrint){print("Calculating confidence intervals can take very long time")}
	        ParamSummary<-.Params.summary(modelParams,EvolModel,EstimationParams$designToEstim,data,1,tmpModelParams$LogLik,nrow(dfData),length(vParFound),RSS,lPrecalculates=lPrecalculates,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik)
    	    }else{ParamSummary<-.Params.summary(modelParams,EvolModel,EstimationParams$designToEstim,data,1,tmpModelParams$LogLik,nrow(dfData),length(vParFound),RSS,lPrecalculates=list(tree.height=lPrecalculates$tree.height))}
	    modelParams<-.cleanUpModelParams(tmpModelParams$modelParams)
	    MaxLikEstim<-list("MaxLik"=vParFound,"ParamSummary"=ParamSummary,"LogLik"=tmpModelParams$LogLik,"ModelParams"=modelParams,"Method.output"=MaximEstim)
	}
    }
    
    if (method=="optim"){## setup likelihood optim
	if(is.null(EstimationParams$StartPoint)){print("No starting position for optimization. Cannot continue")}
	else{
	    par0<-EstimationParams$StartPoint	
	    MaxLikEstim<-optim(
		par=par0,
		fn=function(par){(-1)*.calcLogLikParametrized(par,data,vY,PhylTree,EvolModel,lPrecalculates,EstimationParams,modelParams,TRUE,tol,maxIter,bShouldPrint,minLogLik=minLogLik,parNames=names(par0))}
	    )	    
	}
    }
    if(method=="nlminb"){
	if(is.null(EstimationParams$StartPoint)){print("No starting position for optimization. Cannot continue")}
	else{
	    par0<-EstimationParams$StartPoint
            MaxLikEstim<-nlminb(
                start=par0,
                objective=function(par){(-1)*.calcLogLikParametrized(par,data,vY,PhylTree,EvolModel,lPrecalculates,EstimationParams,modelParams,TRUE,tol,maxIter,bShouldPrint,minLogLik=minLogLik,parNames=names(par0))}
            )
	}
    }
    MaxLikEstim
}

.createCovariancematrix<-function(covmat,n,m,matrixtype,whichcov){
## function creates the covariance matrix from the data provided by the user according to the provided matrix type
## code modified from GLSME.R
## T.F. Hansen, K. Bartoszek; Interpreting the evolutionary regression: the interplay between observational and biological errors in phylogenetic comparative studies; Syst. Biol., 61(3):413-425, 2012
## the GLSME function is not compatible as here we want the matrix to be of the form
## n (species) blocks of size m (traits) while in GLSME m (traits) blocks of size n (species)
## will need a wrapper function if GLSME is included in the mvSLOUCH code for bias correction

    mCov<-matrix(NA,ncol=n*m,nrow=n*m)
    if (is.null(matrixtype)){
        matrixOK<-FALSE
        if ((is.vector(covmat))&&(length(covmat)==1)){matrixtype<-"SingleValue";matrixOK<-TRUE}
        if ((is.vector(covmat))&&(length(covmat)==m)&&(m>1)){matrixtype<-"Vector";matrixOK<-TRUE;}
        if ((is.vector(covmat))&&(length(covmat)==n*m)){matrixtype<-"Diagonal";matrixOK<-TRUE;}
        if ((is.list(covmat))&&(length(covmat)==m)&&(m>1)){matrixtype<-"MatrixList";matrixOK<-TRUE}
        if ((is.matrix(covmat))&&(nrow(covmat)==n*m)&&(ncol(covmat)==n*m)&&(n*m>1)){matrixtype<-"Matrix";matrixOK<-TRUE}
        if (!matrixOK){print(paste("Could not setup ",whichcov," matrix. Please check input or provide a full matrix. Not using a measurement error matrix",sep=""));covmat<-0;matrixtype<-"SingleValue"} ##stop()
    }
    if ((matrixtype=="Matrix")&&(is.matrix(covmat))&&(ncol(covmat)==n*m)&&(nrow(covmat)==n*m)){mCov<-covmat}
    else{
        if (matrixtype=="VectorList"){matrixtype<-"MatrixList"}
        if (matrixtype=="SingleValue"){mCov[1:(n*m),1:(n*m)]<-0;diag(mCov)<-covmat}
        if (matrixtype=="Vector"){mCov[1:(n*m),1:(n*m)]<-0;covmat[which(covmat=="F")]<-0;covmat<-as.numeric(covmat);diagvect<-rep(covmat,n);diag(mCov)<-diagvect}
        if (matrixtype=="Diagonal"){mCov[1:(n*m),1:(n*m)]<-0;diag(mCov)<-covmat}
        if (matrixtype=="CorrelatedPredictors"){mCov<-diag(1,n,n)%x%covmat}
        if (matrixtype=="MatrixList"){
            mCov[1:(n*m),1:(n*m)]<-0;
            for (i in 1:m){
                if ((is.character(covmat[[i]]))&&(covmat[[i]]=="F")){covmat[[i]]<-0};
                if (is.vector(covmat[[i]])){
                    if((length(covmat[[i]])==n)||(length(covmat[[i]])==1)){
                        vcov<-covmat[[i]]
                        covmat[[i]]<-matrix(0,n,n)
                        diag(covmat[[i]])<-vcov
                    }else{print(paste("Could not setup covariance of predictor ",i," of ",whichcov," matrix. Please check input.",sep=""));stop()}
                }
                vindexes<-which((1:(n*m))%%m==i%%m)
                mCov[vindexes,vindexes]<-covmat[[i]]
            }
        }
    }    
    mCov[which(is.na(mCov))]<-0
    if (length(which(is.na(mCov)))>0){print("WARNING: NAs in measurement error covariance matrix changed to 0")}
    mCov<- (mCov+t(mCov))/2    
    EigVals<-eigen(mCov)$values
    if (length(which(EigVals<0)>0)){print(paste("Warning : ",whichcov," matrix has negative eigenvalues!!",sep=""))}
    if (length(which(is.complex(EigVals))>0)){print(paste("Warning : ",whichcov," matrix has complex eigenvalues!!",sep=""))}
    mCov
}
