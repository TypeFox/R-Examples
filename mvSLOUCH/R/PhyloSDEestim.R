.PhyloSDEestim<-function(phyltree,data,kY,regimes=NULL,regimes.times=NULL,root.regime=NULL,predictors=NULL,params=NULL,M.error=NULL){
## predictors have to be given as column numbers
    if (is.null(params)){params<-list()}
    if (!is.element("EvolModel",names(params))){params$EvolModel<-"mvslouch"}
    if (params$EvolModel=="slouch"){params$EvolModel<-"mvslouch";print("Call to slouch not yet implemented, using mvslouch model. You might need to run again with correct input structure.")}	

    phyltree@nodelabels[which(phyltree@nodelabels=="")]<-1:length(which(phyltree@nodelabels==""))
    if(!is.data.frame(data)){
	data<-as.data.frame(data)
	rownames(data)<-phyltree@nodelabels[phyltree@term]
    }
    if (is.null(regimes.times)){regimes.times<-sapply(phyltree@epochs,rev,simplify=FALSE)}
    nullRegimes<-FALSE
    if (is.null(regimes)){
	nullRegimes<-TRUE
	regimes<-sapply(regimes.times,function(reg){rep("reg.1",length(reg)-1)},simplify=FALSE)
    }
    if (!is.list(regimes)){## The regimes are given as a vector in the ouch tree format --- Change them to a list
	if (is.null(root.regime)){root.regime<-regimes[1]}
	vregimes<-regimes
	regimes<-sapply(phyltree@lineages[phyltree@term],function(epch,vregimes){
	    epch<-rev(epch)
	    as.character(sapply(epch[-1],function(reg,vregimes){vregimes[reg]},vregimes=vregimes,simplify=TRUE))
	},vregimes=vregimes,simplify=TRUE)
    }

    bOKregimes<-TRUE
    if ((length(regimes)!=length(regimes.times))||(length(regimes)!=length(phyltree@epochs))){
	bOKregimes<-FALSE;print("Wrong number of regimes, regimes.times vectors")
    }
    else{
	vOKtimes<-sapply(1:length(regimes),function(i,regimes,epochs,regimes.times){
	    bret<-TRUE
	    if (length(regimes[[i]])!=(length(regimes.times[[i]])-1)){bret<-FALSE;print(paste("Problem with ",i,"th regimes or regimes.times",sep=""))}
	    if (length(regimes.times[[i]])<length(epochs[[i]])){bret<-FALSE;print(paste(i,"th regimes.times does not agree with phylogenetic tree",sep=""))}
	    if (length(regimes[[i]])<(length(epochs[[i]])-1)){bret<-FALSE;print(paste(i,"th regimes does not agree with phylogenetic tree",sep=""))}
	    bret	    
	},regimes=regimes,epochs=phyltree@epochs,regimes.times=regimes.times)
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
	if (params$EvolModel=="bm"){kY<-0}
	if (params$EvolModel=="slouch"){kY<-1}
	if ((params$EvolModel=="ouch")&&(is.null(kY)||is.na(kY)||(kY!=ncol(data)))){
	    ##print(paste("Wrong or missing value of kY, ",kY," supplied, changing it to : ",ncol(data),sep=""));
	    kY<-ncol(data)
	}
    	kX<-ncol(data)-kY

	if (!is.element("method",names(params))){params$method<-"glsgc"}
	if (params$EvolModel=="bm"){params$method<-"maxlik"}
	if (!is.element("bShouldPrint",names(params))){params$bShouldPrint<-FALSE}
	if (!is.element("tol",names(params))){params$tol<-.set.tol(params$method)}
	if (!is.element("maxIter",names(params))){params$maxIter<-.set.maxiter(params$method)}
	if (!is.element("maxTries",names(params))){params$maxTries<-10}
	if (!is.element("minLogLik",names(params))){params$minLogLik<- -Inf}
	params$EstimParams<-.set.estimparams(params,kY,kX,length(regimes.types))
	if (params$EstimParams$designToEstim$y0AncState){
	    if (!is.null(root.regime)){params$EstimParams$designToEstim$y0Regime<-which(regimes.types.orig==root.regime)}
	    else{
	     if (is.element("y0Regime",names(params$EstimParams$designToEstim))){params$EstimParams$designToEstim$y0Regime<-which(regimes.types.orig==params$EstimParams$designToEstim$y0Regime)}	
	     else{params$EstimParams$designToEstim$y0Regime<-regimes[[1]][1]} ## need to choose something anyway ...	
	    }
	}

	if (!is.null(predictors)){
	    predictors<-intersect(1:ncol(data),predictors)## if user provided column out of range
	    params$EstimParams$predictors<-predictors
	    predictors<-colnames(data)[predictors]
	}
	if (!is.element("TerminalLabels",params$EstimParams)){params$EstimParams$TerminalLabels<-phyltree@nodelabels[phyltree@term]}
	if (is.null(rownames(data))){rownames(data)<-params$EstimParams$TerminalLabels}
	else{
	    datacols<-colnames(data)
	    data<-data[params$EstimParams$TerminalLabels,]
	    ## the following is needed if we have one variable previous operation will change dfData to a vector
    	     if(!is.data.frame(data)){ data<-as.data.frame(data) ;   rownames(data)<-phyltree@nodelabels[phyltree@term]  ;colnames(data)<-datacols }      	
	} ## order the data correctly according to the tree
	if (params$bShouldPrint){print(Sys.time())}
	lEstimResults<-.maxlik.estim(dfData=data,dfData.Merror=M.error,EvolModel=params$EvolModel,method=params$method,regimes=regimes,regimeTypes=regimes.types,EstimationParams=params$EstimParams,PhylTree=phyltree,regimeTimes=regimes.times,bShouldPrint=params$bShouldPrint,tol=params$tol,maxIter=params$maxIter,maxTries=params$maxTries,minLogLik=params$minLogLik)
	if (params$bShouldPrint){print(Sys.time())}
	lEstimResults<-.correct.names(lEstimResults,regimes.types.orig,if(params$EvolModel!="bm"){colnames(data)[1:kY]}else{NULL},if ((params$EvolModel=="mvslouch")||(params$EvolModel=="bm")||(params$EvolModel=="slouch")){colnames(data)[(kY+1):(kY+kX)]}else{NULL},predictors,params$EvolModel)
	if (!nullRegimes && params$bShouldPrint){print("If you just want the results ignore this message else ... WARNING : if you want to use the returned results in further programs the order of the columns in mPsi variables has to be appropriate to the names of regimes. See the manual for further information.")}
	lEstimResults
    }else{print("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform estimation. Please see manual on their format.")}
}

.correct.names<-function(listobj,vregime.names,vY.names,vX.names,predictors=NULL,EvolModel="mvslouch"){
    if ((EvolModel=="ouch")||(EvolModel=="bm")){
	if (is.element("BrownResult",names(listobj))){listobj$BrownResult<-NULL}
	if (is.element("B",names(listobj))){listobj$B<-NULL}
        if (is.element("B.confidence.interval",names(listobj))){listobj$B.confidence.interval<-NULL}
        if (is.element("B.regression.confidence.interval",names(listobj))){listobj$B.regression.confidence.interval<-NULL}
	if (is.element("Syx",names(listobj))){listobj$Syx<-NULL}
	if (is.element("Sxy",names(listobj))){listobj$Sxy<-NULL}
	if (is.element("Syx.confidence.interval",names(listobj))){listobj$Syx.confidence.interval<-NULL}
	if (is.element("Sxy.confidence.interval",names(listobj))){listobj$Sxy.confidence.interval<-NULL}
	if (is.element("optimal.regression",names(listobj))){listobj$optimal.regression<-NULL}
	if (is.element("A1B",names(listobj))){listobj$A1B<-NULL}
	if (is.element("conditional.cov.matrix",names(listobj))){listobj$conditional.cov.matrix<-NULL}
	if (is.element("conditional.corr.matrix",names(listobj))){listobj$conditional.corr.matrix<-NULL}
	if (is.element("evolutionary.regression",names(listobj))){listobj$evolutionary.regression<-NULL}
	if (is.element("cov.with.optima",names(listobj))){if (length(vY.names)==1){listobj$cov.with.optima<-matrix(listobj$cov.with.optima,1,1)};colnames(listobj$cov.with.optima)<-vY.names;rownames(listobj$cov.with.optima)<-vY.names}
	if (is.element("corr.with.optima",names(listobj))){if (length(vY.names)==1){listobj$corr.with.optima<-matrix(listobj$corr.with.optima,1,1)};colnames(listobj$corr.with.optima)<-vY.names;rownames(listobj$corr.with.optima)<-vY.names}
	if (EvolModel=="bm"){
    	    if (is.element("A",names(listobj))){listobj$A<-NULL}
    	    if (is.element("A.confidence.interval",names(listobj))){listobj$A.confidence.interval<-NULL}
	    if (is.element("mPsi",names(listobj))){listobj$mPsi<-NULL}	
	    if (is.element("mPsi0",names(listobj))){listobj$mPsi0<-NULL}	
	    if (is.element("mPsi.rotated",names(listobj))){listobj$mPsi.rotated<-NULL}	
	    if (is.element("mPsi0.rotated",names(listobj))){listobj$mPsi0.rotated<-NULL}	
    	    if (is.element("expmtA",names(listobj))){listobj$A<-NULL}
	    if (is.element("mPsi.confidence.interval",names(listobj))){listobj$mPsi.confidence.interval<-NULL}	
	    if (is.element("mPsi0.confidence.interval",names(listobj))){listobj$mPsi0.confidence.interval<-NULL}	
	    if (is.element("mPsi.regression.confidence.interval",names(listobj))){listobj$mPsi.regression.confidence.interval<-NULL}	
	    if (is.element("mPsi0.regression.confidence.interval",names(listobj))){listobj$mPsi0.regression.confidence.interval<-NULL}	
    	    if (is.element("expmtA.confidence.interval",names(listobj))){listobj$A.confidence.interval<-NULL}
	}
	if (EvolModel=="ouch"){
	    if (is.element("vX0",names(listobj))){listobj$vX0<-NULL}
	    if (is.element("vX0.confidence.interval",names(listobj))){listobj$vX0.confidence.interval<-NULL}
	    if (is.element("vX0.regression.confidence.interval",names(listobj))){listobj$vX0.regression.confidence.interval<-NULL}
	    if (is.element("Sxx",names(listobj))){listobj$Sxx<-NULL}
	    if (is.element("Sxx.confidence.interval",names(listobj))){listobj$Sxx.confidence.interval<-NULL}
	    if (is.element("optima.correlations",names(listobj))){if (length(vY.names)==1){listobj$optima.correlations<-matrix(listobj$optima.correlations,1,1)};colnames(listobj$optima.correlations)<-vY.names;rownames(listobj$optima.correlations)<-vY.names}    
	}
	if (EvolModel=="bm"){
	    if (is.element("vY0",names(listobj))){listobj$vY0<-NULL}
	    if (is.element("vY0.confidence.interval",names(listobj))){listobj$vY0.confidence.interval<-NULL}
	    if (is.element("vY0.regression.confidence.interval",names(listobj))){listobj$vY0.regression.confidence.interval<-NULL}
	    if (is.element("Syy",names(listobj))){listobj$Syy<-NULL}
	    if (is.element("Syy.confidence.interval",names(listobj))){listobj$Syy.confidence.interval<-NULL}
	}
    }
    if (is.element("A",names(listobj))){if (length(vY.names)==1){listobj$A<-matrix(listobj$A,1,1)};colnames(listobj$A)<-vY.names;rownames(listobj$A)<-vY.names}    
    if (is.element("A.confidence.interval",names(listobj))){colnames(listobj$A.confidence.interval$Lower.end)<-vY.names;rownames(listobj$A.confidence.interval$Lower.end)<-vY.names;colnames(listobj$A.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$A.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$A.confidence.interval$Upper.end)<-vY.names;rownames(listobj$A.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("B",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$B<-matrix(listobj$B,nrow=length(vY.names),ncol=length(vX.names))};colnames(listobj$B)<-vX.names;rownames(listobj$B)<-vY.names}
    if (is.element("B.confidence.interval",names(listobj))){colnames(listobj$B.confidence.interval$Lower.end)<-vX.names;rownames(listobj$B.confidence.interval$Lower.end)<-vY.names;colnames(listobj$B.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$B.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$B.confidence.interval$Upper.end)<-vX.names;rownames(listobj$B.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("B.regression.confidence.interval",names(listobj))){
	if (length(vX.names)==1){rownames(listobj$B.regression.confidence.interval)<-vY.names}    
	else{colnames(listobj$B.regression.confidence.interval$Lower.end)<-vX.names;rownames(listobj$B.regression.confidence.interval$Lower.end)<-vY.names;colnames(listobj$B.regression.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$B.regression.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$B.regression.confidence.interval$Upper.end)<-vX.names;rownames(listobj$B.regression.confidence.interval$Upper.end)<-vY.names}
    }	
    if (is.element("mPsi",names(listobj))){if ((length(vY.names)==1)||(length(vregime.names)==1)){listobj$mPsi<-matrix(listobj$mPsi,nrow=length(vY.names),ncol=length(vregime.names))};colnames(listobj$mPsi)<-vregime.names;rownames(listobj$mPsi)<-vY.names}
    if (is.element("mPsi.confidence.interval",names(listobj))){colnames(listobj$mPsi.confidence.interval$Lower.end)<-vregime.names;rownames(listobj$mPsi.confidence.interval$Lower.end)<-vY.names;colnames(listobj$mPsi.confidence.interval$Estimated.Point)<-vregime.names;rownames(listobj$mPsi.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$mPsi.confidence.interval$Upper.end)<-vregime.names;rownames(listobj$mPsi.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("mPsi.regression.confidence.interval",names(listobj))){
	if (length(vregime.names)==1){rownames(listobj$mPsi.regression.confidence.interval)<-vY.names}    
	else{colnames(listobj$mPsi.regression.confidence.interval$Lower.end)<-vregime.names;rownames(listobj$mPsi.regression.confidence.interval$Lower.end)<-vY.names;colnames(listobj$mPsi.regression.confidence.interval$Estimated.Point)<-vregime.names;rownames(listobj$mPsi.regression.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$mPsi.regression.confidence.interval$Upper.end)<-vregime.names;rownames(listobj$mPsi.regression.confidence.interval$Upper.end)<-vY.names}
    }
    if (is.element("mPsi0",names(listobj))){if (length(vY.names)==1){listobj$mPsi0<-matrix(listobj$mPsi0,nrow=1,ncol=1)};rownames(listobj$mPsi0)<-vY.names}
    if (is.element("mPsi0.confidence.interval",names(listobj))){rownames(listobj$mPsi0.confidence.interval$Lower.end)<-vY.names;rownames(listobj$mPsi0.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$mPsi0.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("mPsi0.regression.confidence.interval",names(listobj))){rownames(listobj$mPsi0.regression.confidence.interval)<-vY.names}   
    if (is.element("vY0",names(listobj))){if (length(vY.names)==1){listobj$vY0<-matrix(listobj$vY0,1,1)};rownames(listobj$vY0)<-vY.names}
    if (is.element("vY0.confidence.interval",names(listobj))){rownames(listobj$vY0.confidence.interval$Lower.end)<-vY.names;rownames(listobj$vY0.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$vY0.confidence.interval$Upper.end)<-vY.names}
    if (is.element("vY0.regression.confidence.interval",names(listobj))){rownames(listobj$vY0.regression.confidence.interval)<-vY.names}
    if (is.element("vX0",names(listobj))){if (length(vX.names)==1){listobj$vX0<-matrix(listobj$vX0,1,1)};rownames(listobj$vX0)<-vX.names}
    if (is.element("vX0.confidence.interval",names(listobj))){rownames(listobj$vX0.confidence.interval$Lower.end)<-vX.names;rownames(listobj$vX0.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$vX0.confidence.interval$Upper.end)<-vX.names}
    if (is.element("vX0.regression.confidence.interval",names(listobj))){rownames(listobj$vX0.regression.confidence.interval)<-vX.names}
    if (is.element("BX0.regression.confidence.interval",names(listobj))){rownames(listobj$BX0.regression.confidence.interval)<-vX.names}
    if (is.element("Syy",names(listobj))){if (length(vY.names)==1){listobj$Syy<-matrix(listobj$Syy,1,1)};colnames(listobj$Syy)<-vY.names;rownames(listobj$Syy)<-vY.names}
    if (is.element("Syy.confidence.interval",names(listobj))){colnames(listobj$Syy.confidence.interval$Lower.end)<-vY.names;rownames(listobj$Syy.confidence.interval$Lower.end)<-vY.names;colnames(listobj$Syy.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$Syy.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$Syy.confidence.interval$Upper.end)<-vY.names;rownames(listobj$Syy.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("Syx",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$Syx<-matrix(listobj$Syx,nrow=length(vY.names),ncol=length(vX.names))};colnames(listobj$Syx)<-vX.names;rownames(listobj$Syx)<-vY.names}
    if (is.element("Syx.confidence.interval",names(listobj))){colnames(listobj$Syx.confidence.interval$Lower.end)<-vX.names;rownames(listobj$Syx.confidence.interval$Lower.end)<-vY.names;colnames(listobj$Syx.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$Syx.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$Syx.confidence.interval$Upper.end)<-vX.names;rownames(listobj$Syx.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("Sxy",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$Sxy<-matrix(listobj$Sxy,ncol=length(vY.names),nrow=length(vX.names))};colnames(listobj$Sxy)<-vY.names;rownames(listobj$Sxy)<-vX.names}
    if (is.element("Sxy.confidence.interval",names(listobj))){colnames(listobj$Sxy.confidence.interval$Lower.end)<-vY.names;rownames(listobj$Sxy.confidence.interval$Lower.end)<-vX.names;colnames(listobj$Sxy.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$Sxy.confidence.interval$Estimated.Point)<-vX.names;colnames(listobj$Sxy.confidence.interval$Upper.end)<-vY.names;rownames(listobj$Sxy.confidence.interval$Upper.end)<-vX.names}    
    if (is.element("Sxx",names(listobj))){if (length(vX.names)==1){listobj$Sxx<-matrix(listobj$Sxx,1,1)};colnames(listobj$Sxx)<-vX.names;rownames(listobj$Sxx)<-vX.names}
    if (is.element("Sxx.confidence.interval",names(listobj))){colnames(listobj$Sxx.confidence.interval$Lower.end)<-vX.names;rownames(listobj$Sxx.confidence.interval$Lower.end)<-vX.names;colnames(listobj$Sxx.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$Sxx.confidence.interval$Estimated.Point)<-vX.names;colnames(listobj$Sxx.confidence.interval$Upper.end)<-vX.names;rownames(listobj$Sxx.confidence.interval$Upper.end)<-vX.names}    
    if (is.element("expmtA",names(listobj))){if (length(vY.names)==1){listobj$expmtA<-matrix(listobj$expmtA,1,1)};colnames(listobj$expmtA)<-vY.names;rownames(listobj$expmtA)<-vY.names}
    if (is.element("mPsi.rotated",names(listobj))){if ((length(vY.names)==1)||(length(vregime.names)==1)){listobj$mPsi.rotated<-matrix(listobj$mPsi.rotated,nrow=length(vY.names),ncol=length(vregime.names))};colnames(listobj$mPsi.rotated)<-vregime.names;rownames(listobj$mPsi.rotated)<-vY.names}
    if (is.element("mPsi0.rotated",names(listobj))){if (length(vY.names)==1)   {listobj$mPsi0.rotated<-matrix(listobj$mPsi0.rotated,nrow=1,ncol=1)};rownames(listobj$mPsi0.rotated)<-vY.names}
    if (is.element("optimal.regression",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$optimal.regression<-matrix(listobj$optimal.regression,nrow=length(vY.names),ncol=length(vX.names))};rownames(listobj$optimal.regression)<-vY.names;colnames(listobj$optimal.regression)<-vX.names}
    if (is.element("A1B",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$A1B<-matrix(listobj$A1B,nrow=length(vY.names),ncol=length(vX.names))};rownames(listobj$A1B)<-vY.names;colnames(listobj$A1B)<-vX.names}
    if (is.element("cov.matrix",names(listobj))){colnames(listobj$cov.matrix)<-c(vY.names,vX.names);rownames(listobj$cov.matrix)<-c(vY.names,vX.names)}
    if (is.element("corr.matrix",names(listobj))){colnames(listobj$corr.matrix)<-c(vY.names,vX.names);rownames(listobj$corr.matrix)<-c(vY.names,vX.names)}
    if (is.element("conditional.cov.matrix",names(listobj))){if (length(vY.names)==1){listobj$conditional.cov.matrix<-matrix(listobj$conditional.cov.matrix,1,1)};colnames(listobj$conditional.cov.matrix)<-vY.names;rownames(listobj$conditional.cov.matrix)<-vY.names}
    if (is.element("conditional.corr.matrix",names(listobj))){if (length(vY.names)==1){listobj$conditional.corr.matrix<-matrix(listobj$conditional.corr.matrix,1,1)};colnames(listobj$conditional.corr.matrix)<-vY.names;rownames(listobj$conditional.corr.matrix)<-vY.names}
    if (is.element("evolutionary.regression",names(listobj))){
	if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$evolutionary.regression<-matrix(listobj$evolutionary.regression,ncol=length(vX.names),nrow=length(vY.names))}
	if (is.null(predictors)){colnames(listobj$evolutionary.regression)<-vX.names;rownames(listobj$evolutionary.regression)<-vY.names}
	else{
	    if (length(setdiff(predictors,vX.names))>0){
		listobj$regression.Y.on.X<-listobj$evolutionary.regression
		colnames(listobj$regression.Y.on.X)<-vX.names;rownames(listobj$regression.Y.on.X)<-vY.names
		listobj$evolutionary.regression<-NULL
		tryCatch({
		    listobj$evolutionary.regression<-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),predictors]%*%solve(listobj$cov.matrix[predictors,predictors])		    		
		    if ((length(predictors)==1)||(length(setdiff(c(vY.names,vX.names),predictors))==1)){listobj$evolutionary.regression<-matrix(listobj$evolutionary.regression,ncol=length(predictors),nrow=length(setdiff(c(vY.names,vX.names),predictors)))}
		    colnames(listobj$evolutionary.regression)<-predictors;rownames(listobj$evolutionary.regression)<-setdiff(c(vY.names,vX.names),predictors)
		},error=function(e){print(paste("Error in evolutionary regression calculation",e))})		
		if (is.element("conditional.cov.matrix",names(listobj))){
		    listobj$conditional.cov.Y.on.X<-listobj$conditional.cov.matrix
		    listobj$conditional.corr.Y.on.X<-listobj$conditional.corr.matrix
		}
	       listobj$conditional.cov.matrix<-NULL
	       listobj$conditional.corr.matrix<-NULL
		tryCatch({
		    listobj$conditional.cov.matrix<-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),predictors]%*%solve(listobj$cov.matrix[predictors,predictors])%*%listobj$cov.matrix[predictors,setdiff(c(vY.names,vX.names),predictors)]
		    listobj$conditional.corr.matrix<-cov2cor(listobj$conditional.cov.matrix) ##produces matrix, no need to correct if 1-dim
		    if (length(setdiff(c(vY.names,vX.names),predictors))==1){listobj$conditional.cov.matrix<-matrix(listobj$conditional.cov.matrix,ncol=1,nrow=1)}
		    colnames(listobj$conditional.cov.matrix)<-setdiff(c(vY.names,vX.names),predictors);rownames(listobj$conditional.cov.matrix)<-setdiff(c(vY.names,vX.names),predictors)
		    colnames(listobj$conditional.corr.matrix)<-setdiff(c(vY.names,vX.names),predictors);rownames(listobj$conditional.corr.matrix)<-setdiff(c(vY.names,vX.names),predictors)
		},error=function(e){print(paste("Error in conditional covariance calculation",e))})		
	    }	    
	    else{colnames(listobj$evolutionary.regression)<-vX.names;rownames(listobj$evolutionary.regression)<-vY.names}	
	}
    }else{
	if ((is.element("cov.matrix",names(listobj)))&&(!is.null(predictors))){
	    listobj$evolutionary.regression<-NULL
	    tryCatch({	
    		listobj$evolutionary.regression<-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),predictors]%*%solve(listobj$cov.matrix[predictors,predictors])
		    if ((length(predictors)==1)||(length(setdiff(c(vY.names,vX.names),predictors))==1)){listobj$evolutionary.regression<-matrix(listobj$evolutionary.regression,ncol=length(predictors),nrow=length(setdiff(c(vY.names,vX.names),predictors)))}
		    colnames(listobj$evolutionary.regression)<-predictors;rownames(listobj$evolutionary.regression)<-setdiff(c(vY.names,vX.names),predictors)
	    },error=function(e){print(paste("Error in evolutionary regression calculation",e))})		
	    if (is.element("conditional.cov.matrix",names(listobj))){
		listobj$conditional.cov.Y.on.X<-listobj$conditional.cov.matrix
		listobj$conditional.corr.Y.on.X<-listobj$conditional.corr.matrix
	    }
	    listobj$conditional.cov.matrix<-NULL
	    listobj$conditional.corr.matrix<-NULL
	    tryCatch({
		listobj$conditional.cov.matrix<-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),predictors]%*%solve(listobj$cov.matrix[predictors,predictors])%*%listobj$cov.matrix[predictors,setdiff(c(vY.names,vX.names),predictors)]
		listobj$conditional.corr.matrix<-cov2cor(listobj$conditional.cov.matrix) ##produces matrix, no need to correct if 1-dim
		if (length(setdiff(c(vY.names,vX.names),predictors))==1){listobj$conditional.cov.matrix<-matrix(listobj$conditional.cov.matrix,ncol=1,nrow=1)}
		colnames(listobj$conditional.cov.matrix)<-setdiff(c(vY.names,vX.names),predictors);rownames(listobj$conditional.cov.matrix)<-setdiff(c(vY.names,vX.names),predictors)
		colnames(listobj$conditional.corr.matrix)<-setdiff(c(vY.names,vX.names),predictors);rownames(listobj$conditional.corr.matrix)<-setdiff(c(vY.names,vX.names),predictors)
	    },error=function(e){print(paste("Error in conditional covariance calculation",e))})		
	}
    }
    if (is.element("trait.regression",names(listobj))){
        listobj$trait.regression<-sapply(1:length(vY.names),function(i,vY.names,vX.names,trait.reg){
	    		    trait.reg<-matrix(trait.reg[[i]],nrow=1);colnames(trait.reg)<-c(vY.names[-i],vX.names);rownames(trait.reg)<-vY.names[i]
			    trait.reg
	},vY.names=vY.names,vX.names=vX.names,trait.reg=listobj$trait.regression,simplify=FALSE)
        if (length(setdiff(predictors,vX.names))>0){
	    listobj$Y.regression<-listobj$trait.regression
	    if (is.element("cov.matrix",names(listobj))){
		vtraits<-setdiff(c(vY.names,vX.names),predictors)
		vtraitsIndex<-sort(sapply(vtraits,function(trait,vNames){which(vNames==trait)},vNames=c(vY.names,vX.names),simplify=TRUE))
		listobj$trait.regression<-NULL
		tryCatch({		
		    listobj$trait.regression<-sapply(1:length(vtraits),function(i,vtraitsIndex,mCov){
				mCov[vtraitsIndex[i],-vtraitsIndex[i]]%*%solve(mCov[-vtraitsIndex[i],-vtraitsIndex[i]])    
		    },vtraitsIndex=vtraitsIndex,mCov=listobj$cov.matrix,simplify=FALSE)
		    listobj$trait.regression<-sapply(1:length(vtraits),function(i,vY.names,vX.names,vtraitsIndex,trait.reg){
			    trait.reg<-matrix(trait.reg[[i]],nrow=1);colnames(trait.reg)<-c(vY.names,vX.names)[-vtraitsIndex[i]];rownames(trait.reg)<-c(vY.names,vX.names)[vtraitsIndex[i]]
			    trait.reg
		    },vY.names=vY.names,vX.names=vX.names,vtraitsIndex=vtraitsIndex,trait.reg=listobj$trait.regression,simplify=FALSE)
		},error=function(e){print(paste("Error in trait regression calculation",e))})		
	    }else{listobj$trait.regression<-NULL}
	}
    }
    if (is.element("cov.with.optima",names(listobj))){if (length(vY.names)==1){listobj$cov.with.optima<-matrix(listobj$cov.with.optima,1,1)};colnames(listobj$cov.with.optima)<-vY.names;rownames(listobj$cov.with.optima)<-vY.names}
    if (is.element("corr.with.optima",names(listobj))){if (length(vY.names)==1){listobj$corr.with.optima<-matrix(listobj$corr.with.optima,1,1)};colnames(listobj$corr.with.optima)<-vY.names;rownames(listobj$corr.with.optima)<-vY.names}
    if (is.element("optima.cov.matrix",names(listobj))){if (length(vY.names)==1){listobj$optima.cov.matrix<-matrix(listobj$optima.cov.matrix,1,1)};colnames(listobj$optima.cov.matrix)<-vY.names;rownames(listobj$optima.cov.matrix)<-vY.names}
    if (is.element("optima.corr.matrix",names(listobj))){if (length(vY.names)==1){listobj$optima.corr.matrix<-matrix(listobj$optima.corr.matrix,1,1)};colnames(listobj$optima.corr.matrix)<-vY.names;rownames(listobj$optima.corr.matrix)<-vY.names}
    if (is.element("stationary.cov.matrix",names(listobj))){
	if (length(vY.names)==1){listobj$stationary.cov.matrix<-matrix(listobj$stationary.cov.matrix,1,1)};colnames(listobj$stationary.cov.matrix)<-vY.names;rownames(listobj$stationary.cov.matrix)<-vY.names
	if (length(vY.names)==1){listobj$stationary.corr.matrix<-matrix(listobj$stationary.corr.matrix,1,1)};colnames(listobj$stationary.corr.matrix)<-vY.names;rownames(listobj$stationary.corr.matrix)<-vY.names
	if ((length(vX.names)==0)&&(!is.null(predictors))){
	    vResps<-setdiff(vY.names,predictors)
	    tryCatch({
		if (!is.element("limiting.regression",names(listobj))){
		    listobj$limiting.regression<-listobj$stationary.cov.matrix[vResps,predictors]%*%solve(listobj$stationary.cov.matrix[predictors,predictors])
		    if ((length(vResps)==1)||(length(predictors)==1)){listobj$limiting.regression<-matrix(listobj$limiting.regression,nrow=length(vResps),ncol=length(predictors))};rownames(listobj$limiting.regression)<-vResps;colnames(listobj$limiting.regression)<-predictors
		}
		if (!is.element("limiting.relationship.cov.matrix",names(listobj))){
		    listobj$limiting.relationship.cov.matrix<-listobj$limiting.regression%*%listobj$stationary.cov.matrix[predictors,predictors]%*%t(listobj$limiting.regression)
		    listobj$limiting.relationship.corr.matrix<-cov2cor(listobj$limiting.relationship.cov.matrix)
		    if (length(vResps)==1){listobj$limiting.relationship.cov.matrix<-matrix(listobj$limiting.relationship.cov.matrix,1,1)};colnames(listobj$limiting.relationship.cov.matrix)<-vResps;rownames(listobj$limiting.relationship.cov.matrix)<-vResps    
		}
		if (!is.element("cov.with.limit",names(listobj))){
		    listobj$cov.with.limit<-listobj$stationary.cov.matrix[vResps,predictors]%*%t(listobj$limiting.regression)
		    mcov.curr<-listobj$stationary.cov.matrix[vResps,vResps]
		    mcov.stat<-listobj$limiting.relationship.cov.matrix[vResps,vResps]
		    if (length(vResps)==1){
			mcov.curr<-matrix(mcov.curr,1,1)
			mcov.stat<-matrix(mcov.stat,1,1)
		    }
		    listobj$corr.with.limit<-apply(matrix(0:((length(vResps))^2-1),length(vResps),length(vResps),byrow=TRUE),c(1,2),function(ij,kY,mcov.curr,mcov.stat,mcov.with){i<-ij%/%kY+1;j<-ij%%kY+1;mcov.with[i,j]/(sqrt(mcov.curr[i,i]*mcov.stat[j,j]))},kY=length(vResps),mcov.curr=mcov.curr,mcov.stat=mcov.stat,mcov.with=listobj$cov.with.limit)
		    if (length(vResps)==1){listobj$cov.with.limit<-matrix(listobj$cov.with.limit,1,1)};colnames(listobj$cov.with.limit)<-vResps;rownames(listobj$cov.with.limit)<-vResps    
		}
	    },error=function(e){print(paste("Error in limiting regression and covariance calculation",e))})		
	}
	tryCatch({    
	    if ((!is.element("limiting.trait.regression",names(listobj)))&&(length(vY.names)>1)){ ## no point in doing the regressions if there is only a single trait
		listobj$limiting.trait.regression<-sapply(1:length(vY.names),function(i,mCov,vY.names){
			    limiting.trait.reg<-mCov[i,-i]%*%solve(mCov[-i,-i])    
			    limiting.trait.reg<-matrix(limiting.trait.reg,nrow=1);colnames(limiting.trait.reg)<-vY.names[-i];rownames(limiting.trait.reg)<-vY.names[i]
			    limiting.trait.reg	
		},mCov=listobj$stationary.cov.matrix,vY.names=vY.names,simplify=FALSE)
	    }
	},error=function(e){print(paste("Error in limiting trait regression calculation",e))})
    }   
    if (is.element("StS",names(listobj))){colnames(listobj$StS)<-c(vY.names,vX.names);rownames(listobj$StS)<-c(vY.names,vX.names)}
    if (is.element("lower.summary",names(listobj))){
        if (is.element("LogLik",names(listobj$lower.summary))){listobj$lower.summary$LogLik<-NULL}
	if (is.element("dof",names(listobj$lower.summary))){listobj$lower.summary$dof<-NULL}
	if (is.element("m2loglik",names(listobj$lower.summary))){listobj$lower.summary$m2loglik<-NULL}
	if (is.element("aic",names(listobj$lower.summary))){listobj$lower.summary$aic<-NULL}
	if (is.element("aic.c",names(listobj$lower.summary))){listobj$lower.summary$aic.c<-NULL}
	if (is.element("sic",names(listobj$lower.summary))){listobj$lower.summary$sic<-NULL}
	if (is.element("bic",names(listobj$lower.summary))){listobj$lower.summary$bic<-NULL}
	if (is.element("RSS",names(listobj$lower.summary))){listobj$lower.summary$RSS<-NULL}
    }
    if (is.element("upper.summary",names(listobj))){
        if (is.element("LogLik",names(listobj$upper.summary))){listobj$upper.summary$LogLik<-NULL}
	if (is.element("dof",names(listobj$upper.summary))){listobj$upper.summary$dof<-NULL}
	if (is.element("m2loglik",names(listobj$upper.summary))){listobj$upper.summary$m2loglik<-NULL}
	if (is.element("aic",names(listobj$upper.summary))){listobj$upper.summary$aic<-NULL}
	if (is.element("aic.c",names(listobj$upper.summary))){listobj$upper.summary$aic.c<-NULL}
	if (is.element("sic",names(listobj$upper.summary))){listobj$upper.summary$sic<-NULL}
	if (is.element("bic",names(listobj$upper.summary))){listobj$upper.summary$bic<-NULL}
	if (is.element("RSS",names(listobj$upper.summary))){listobj$upper.summary$RSS<-NULL}
    }
    if (is.element("regressCovar",names(listobj))){listobj$regressCovar<-NULL}    
    sapply(listobj,function(obj,vregime.names,vY.names,vX.names,predictors,EvolModel){if (is.list(obj)){
    .correct.names(obj,vregime.names,vY.names,vX.names,predictors,EvolModel)}else{obj}},vregime.names=vregime.names,vY.names=vY.names,vX.names=vX.names,predictors=predictors,EvolModel=EvolModel,simplify=FALSE)
}

.set.tol<-function(method){
    tol=switch(method,
	glsgc=c(0.01,0.001),
	gridigls=0.001,
	igls=0.001,
	grid=NA,
	gridgls=NA,
	gls=NA,
	maxlik=NA,
	optim=NA,
	nlminb=NA,
	EM=NA,
	MCEM=NA	
    )
    tol
}

.set.maxiter<-function(method){
    maxiter=switch(method,
	glsgc=c(10,50),
	gridigls=50,
	igls=50,
	grid=NA,
	gridgls=NA,
	gls=NA,
	maxlik=NA,
	optim=NA,
	nlminb=NA,
	EM=NA,
	MCEM=NA	
    )
    maxiter
}

.set.estimparams<-function(params,kY,kX,numregs){
    if (!is.element("EstimParams",names(params))){EstimParams<-list()}
    else{EstimParams<-params$EstimParams}
    EstimParams$vVars<-NULL
    EstimParams$conditional<-FALSE
    if (!is.element("kY",names(EstimParams))){EstimParams$kY<-kY}
    if (!is.element("kX",names(EstimParams))){EstimParams$kX<-kX}
    if (!is.element("Atype",names(EstimParams))){EstimParams$Atype<-"Invertible"}
    if (!is.element("diagA",names(EstimParams))){EstimParams$diagA<-NULL}
    if (!is.element("Syytype",names(EstimParams))){EstimParams$Syytype<-"UpperTri"}
    if (!is.element("diagSyy",names(EstimParams))){EstimParams$diagSyy<-"Positive"}
    if (!is.element("signsA",names(EstimParams))){EstimParams$signsA<-NULL}
    if (!is.element("signsB",names(EstimParams))){EstimParams$signsB<-NULL}
    if (!is.element("signsmPsi",names(EstimParams))){EstimParams$signsmPsi<-NULL}
    if (!is.element("signsmPsi0",names(EstimParams))){EstimParams$signsmPsi0<-NULL}
    if (!is.element("signsvY0",names(EstimParams))){EstimParams$vY0<-NULL}
    if (!is.element("signsvX0",names(EstimParams))){EstimParams$vX0<-NULL}
    if (!is.element("signsSyy",names(EstimParams))){EstimParams$Syy<-NULL}    
    if (!is.element("signsSyx",names(EstimParams))){EstimParams$Syx<-NULL}
    if (!is.element("signsSxy",names(EstimParams))){EstimParams$Sxy<-NULL}
    if (!is.element("signsSxx",names(EstimParams))){EstimParams$Sxx<-NULL}        
    if (!is.element("maximMethod",names(EstimParams))){EstimParams$maximMethod<-"optim"}
    if (!is.element("conf.level",names(EstimParams))){EstimParams$conf.level<-0.95}
    if (!is.element("calcCI",names(EstimParams))){EstimParams$calcCI<-FALSE}
    if (!is.element("designToEstim",names(EstimParams))){EstimParams$designToEstim<-list()}
    if (!is.element("y0",names(EstimParams$designToEstim))){EstimParams$designToEstim$y0<-TRUE}
    if (!is.element("psi",names(EstimParams$designToEstim))){EstimParams$designToEstim$psi<-TRUE}
    if (!is.element("psi0",names(EstimParams$designToEstim))){EstimParams$designToEstim$psi0<-FALSE}
    if ((!is.element("X0",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$X0<-FALSE}
    if ((!is.element("B",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$B<-TRUE}
    if ((!is.element("BX0",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$BX0<-FALSE}
    if ((!is.element("UseX0",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$UseX0<-TRUE}	
    if (!is.element("y0AncState",names(EstimParams$designToEstim))||(numregs==1)){EstimParams$designToEstim$y0AncState<-TRUE}
    if ((!is.element("y0OnlyFixed",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$y0OnlyFixed<-FALSE}
    if (!is.element("SimpReg",names(EstimParams$designToEstim))){EstimParams$designToEstim$SimpReg<-FALSE}    	
    if ((!is.element("iRegLin",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$iRegLin<-TRUE}
    if ((!is.element("optimMethod",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$optimMethod<-"optim"}
    if ((!is.element("BFullXNA",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$BFullXNA<-TRUE}	
    if ((!is.element("FullNAYX",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$FullNAYX<-TRUE}	
    if (!is.element("FullNAY",names(EstimParams$designToEstim))){EstimParams$designToEstim$FullNAY<-TRUE}	
    if (!is.element("sigmaRule",names(EstimParams$designToEstim))){EstimParams$designToEstim$sigmaRule<-3}	
    if (!is.element("Fixed",names(EstimParams))){EstimParams$Fixed<-list()}
    if (!is.element("mPsi",names(EstimParams$Fixed))){EstimParams$Fixed$mPsi<-NA}
    if (!is.element("mPsi0",names(EstimParams$Fixed))){	if ((kY>0)&&((numregs==1)||(!EstimParams$designToEstim$psi0)||(!EstimParams$designToEstim$y0AncState))){EstimParams$Fixed$mPsi0<-matrix(0,ncol=1,nrow=kY)}else{EstimParams$Fixed$mPsi0<-NA}}
    ## No need to correct for estimating both mPsi0 and Y0 if there is only one regime as this mPsi0 is immediately zeroed in one regime case
    if (!is.element("vY0",names(EstimParams$Fixed))){EstimParams$Fixed$vY0<-NA}
    if (!is.element("vX0",names(EstimParams$Fixed))){EstimParams$Fixed$vX0<-NA}
    if (!is.element("Sxx",names(EstimParams$Fixed))){EstimParams$Fixed$Sxx<-NA}
    if (!is.element("B",names(EstimParams$Fixed))){EstimParams$Fixed$B<-NA}
    if (!is.element("Sxy",names(EstimParams$Fixed))){if((kX>0)&&(kY>0)){EstimParams$Fixed$Sxy<-matrix(0,nrow=kX,ncol=kY)}else{EstimParams$Fixed$Sxy<-NA}}    
    if (!is.element("Syx",names(EstimParams$Fixed))){if((kX>0)&&(kY>0)){EstimParams$Fixed$Syx<-matrix(0,nrow=kY,ncol=kX)}else{EstimParams$Fixed$Syx<-NA}}        
    if (!is.element("KnownParams",names(EstimParams))){
	vKnownNames<-c()
	EstimParams$KnownParams<-list()
	j<-1	
	#if ((length(EstimParams$Fixed)>0)&&(length(EstimParams$KnownParams)>0)){
	if (length(EstimParams$Fixed)>0){#&&(length(EstimParams$KnownParams)>0)){
	    for (i in 1:length(EstimParams$Fixed)){
		if (!is.na(EstimParams$Fixed[[i]][1])){
		    EstimParams$KnownParams[[j]]<-EstimParams$Fixed[[i]]
		    vKnownNames<-c(vKnownNames,names(EstimParams$Fixed)[i])
		    j<-j+1
		}	    
	    }
	}
	if (j>1){names(EstimParams$KnownParams)<-vKnownNames}
    }else{
	vListNull<-c()
	if ((length(EstimParams$Fixed)>0)&&(length(EstimParams$KnownParams)>0)){	
	    for (j in 1:length(EstimParams$KnownParams)){
		i<-which(names(EstimParams$Fixed)==names(EstimParams$KnownParams)[j])
		if (is.na(EstimParams$Fixed[[i]][1])){vListNull<-c(vListNull,names(EstimParams$KnowParams)[j])}
	    }	
	}
	if (length(vListNull)>0){for(par in vListNull){EstimParams$KnownParams[[which(names(EstimParams$KnownParams)==par)]]<-NULL}}
    }
    if ((!is.element("StartPoint",names(EstimParams)))&&(params$method!="grid")&&(params$method!="gridgls")&&(params$method!="gridigls")){
	parnames<-.set.paramatrizationnames(EstimParams,params$EvolModel,kY,kX,numregs)
	EstimParams$StartPoint<-rnorm(length(parnames))
	names(EstimParams$StartPoint)<-parnames
    }
    if ((!is.element("mGrid",names(EstimParams)))&&((params$method=="grid")||(params$method=="gridgls")||(params$method=="gridigls"))){
	if (!is.element("RegularGrid",names(EstimParams))){EstimParams$RegularGrid<-FALSE}
	if (!EstimParams$RegularGrid){
	    parnames<-.set.paramatrizationnames(EstimParams,params$EvolModel,kY,kX,numregs)
	    EstimParams$mGrid<-matrix(rnorm(length(parnames)*1000),length(parnames))
	    colnames(EstimParams$mGrid)<-parnames
	    EstimParams$mGrid<-.resizegrid(EstimParams)
	}
    }               
    EstimParams
}

.set.paramatrizationnames<-function(EstimParams,EvolModel,kY,kX,numregs){
    parnames<-c()
    if (EvolModel=="bm"){
    	if (!is.element("vY0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("vY0",kY))}
    	if (!is.element("Sxx",names(EstimParams$Fixed))){
	    if (is.element("Sxxtype",names(EstimParams))){
		Sxxlength=switch(EstimParams$Sxxtype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     Symmetric={(kY+1)*kY/2},
                     UpperTri={(kY+1)*kY/2},
                     LowerTri={(kY+1)*kY/2},
                     Any={kY^2}
		)
		parnames<-c(parnames,.generatenames("Sxx",Sxxlength))
	    }
	}	

    }
    if (EvolModel=="ouch"){
    	if (!is.element("A",names(EstimParams$Fixed))){
	    if (is.element("Atype",names(EstimParams))){
		Alength=switch(EstimParams$Atype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     SymmetricPositiveDefinite={(kY+1)*kY/2},
                     Symmetric={(kY+1)*kY/2},
                     TwoByTwo={4},
                     UpperTri={(kY+1)*kY/2},
                     LowerTri={(kY+1)*kY/2},
                     DecomposablePositive={kY^2},
                     DecomposableNegative={kY^2},
                     DecomposableReal={kY^2},
                     Invertible={kY^2}
		)
		parnames<-c(parnames,.generatenames("A",Alength))
	    }
	}
	if (!is.element("vY0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("vY0",kY))}
	if (!is.element("mPsi",names(EstimParams$Fixed))){
	    if (is.element("mPsitype",names(EstimParams))){
		psilength=switch(EstimParams$mPsitype,
	    	     Global={kY},
	    	     Regimes={kY*numregs}
		)
		parnames<-c(parnames,.generatenames("Psi",psilength))
	    }
	}
	if (!is.element("mPsi0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("Psi0",kY))}
	if (!is.element("Syy",names(EstimParams$Fixed))){
	    if (is.element("Syytype",names(EstimParams))){
		Syylength=switch(EstimParams$Syytype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     Symmetric={(kY+1)*kY/2},
                     UpperTri={(kY+1)*kY/2},
                     LowerTri={(kY+1)*kY/2},
                     Any={kY^2}
		)
		parnames<-c(parnames,.generatenames("Syy",Syylength))
	    }
	}	   
    }
    if (EvolModel=="slouch"){}
    if (EvolModel=="mvslouch"){
	if (!is.element("A",names(EstimParams$Fixed))){
	    if (is.element("Atype",names(EstimParams))){
		Alength=switch(EstimParams$Atype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     SymmetricPositiveDefinite={(kY+1)*kY/2},
                     Symmetric={(kY+1)*kY/2},
                     TwoByTwo={4},
                     UpperTri={(kY+1)*kY/2},
                     LowerTri={(kY+1)*kY/2},
                     DecomposablePositive={kY^2},
                     DecomposableNegative={kY^2},
                     DecomposableReal={kY^2},
                     Invertible={kY^2}
		)
		parnames<-c(parnames,.generatenames("A",Alength))
	    }
	}
	if (!is.element("B",names(EstimParams$Fixed))){
	    if (is.element("Btype",names(EstimParams))){
		Blength=switch(EstimParams$Btype,
	    	     MinusA={0},
	    	     SingleValue={1},
	    	     SingleValueDiagonal={1},
                     Diagonal={min(kY,kX)},
                     Symmetric={min(kY,kX)*(min(kY,kX)+1)/2+kY*kX-min(kY,kX)^2},
                     Any={kY*kX}
		)
		parnames<-c(parnames,.generatenames("B",Blength))
	    }
	}
	if (!is.element("vX0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("vX0",kX))}	
	if (!is.element("vY0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("vY0",kY))}
	if (!is.element("mPsi",names(EstimParams$Fixed))){
	    if (is.element("mPsitype",names(EstimParams))){
		psilength=switch(EstimParams$mPsitype,
	    	     Global={kY},
	    	     Regimes={kY*numregs}
		)
		parnames<-c(parnames,.generatenames("Psi",psilength))
	    }
	}
	if (!is.element("mPsi0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("mPsi0",kY))}
	if (!is.element("Syy",names(EstimParams$Fixed))){
	    if (is.element("Syytype",names(EstimParams))){
		Syylength=switch(EstimParams$Syytype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     Symmetric={(kY+1)*kY/2},
                     UpperTri={(kY+1)*kY/2},
                     LowerTri={(kY+1)*kY/2},
                     Any={kY^2}
		)
		parnames<-c(parnames,.generatenames("Syy",Syylength))
	    }
	}	
	if (!is.element("Syx",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("Syx",kY*kX))}	
	if (!is.element("Sxy",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("Sxy",kY*kX))}	
	if (!is.element("Sxx",names(EstimParams$Fixed))){
	    if (is.element("Sxxtype",names(EstimParams))){
		Sxxlength=switch(EstimParams$Sxxtype,
		     OneSigmaDiagonal={1},
                     Diagonal={kX},
                     Symmetric={(kX+1)*kX/2},
                     Any={kX^2}
		)
		parnames<-c(parnames,.generatenames("Sxx",Sxxlength))
	    }
	}	
    }    
    parnames
}

.generatenames<-function(prefix,len){
    gennenames<-c()	
    if (len==1){gennednames<-c(prefix)}
    if (len==2){gennednames<-c(paste(prefix,"start",sep=""),paste(prefix,"end",sep=""))}
    if (len>2){
	gennednames<-sapply(1:len,function(x,prefix){paste(prefix,"_",x,sep="")},prefix=prefix)
	gennednames[1]<-paste(prefix,"start",sep="")
	gennednames[len]<-paste(prefix,"end",sep="")	                             
    }	    	    
    gennednames
}

.resizegrid<-function(EstimParams){
## if a random grid was generated resize it into values between allowed min and max for each parameter
    if (is.element("Aminmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="A")
	if (length(colstr)==0){	colstr<-(which(colnames(EstimParams$mGrid)=="Astart")):(which(colnames(EstimParams$mGrid)=="Aend"))}
	for (i in colstr){
	    minVal<-EstimParams$Aminmax[which(rownames(EstimParams$Aminmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$Aminmax[which(rownames(EstimParams$Aminmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("Bminmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="B")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="Bstart")):(which(colnames(EstimParams$mGrid)=="Bend"))}
	for (i in colstr){
	    minVal<-EstimParams$Bminmax[which(rownames(EstimParams$Bminmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$Bminmax[which(rownames(EstimParams$Bminmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("mPsiminmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="Psi")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="Psistart")):(which(colnames(EstimParams$mGrid)=="Psiend"))}
	for (i in colstr){
	    minVal<-EstimParams$mPsiminmax[which(rownames(EstimParams$mPsiminmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$mPsiminmax[which(rownames(EstimParams$mPsiminmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("mPsi0minmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="mPsi0")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="Psi0start")):(which(colnames(EstimParams$mGrid)=="Psi0end"))}
	for (i in colstr){
	    minVal<-EstimParams$Psi0minmax[which(rownames(EstimParams$Psi0minmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$Psi0minmax[which(rownames(EstimParams$Psi0minmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("vY0minmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="vY0")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="vY0start")):(which(colnames(EstimParams$mGrid)=="vY0end"))}
	for (i in colstr){
	    minVal<-EstimParams$vY0minmax[which(rownames(EstimParams$vY0minmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$vY0minmax[which(rownames(EstimParams$vY0minmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("vX0minmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="vX0")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="vX0start")):(which(colnames(EstimParams$mGrid)=="vX0end"))}
	for (i in colstr){
	    minVal<-EstimParams$vX0minmax[which(rownames(EstimParams$vX0minmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$vX0minmax[which(rownames(EstimParams$vX0minmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("Syyminmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="Syy")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="Syystart")):(which(colnames(EstimParams$mGrid)=="Syyend"))}
	for (i in colstr){
	    minVal<-EstimParams$Syyminmax[which(rownames(EstimParams$Syyminmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$Syyminmax[which(rownames(EstimParams$Syyminmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("Syxminmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="Syx")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="Syxstart")):(which(colnames(EstimParams$mGrid)=="Syxend"))}
	for (i in colstr){
	    minVal<-EstimParams$Syxminmax[which(rownames(EstimParams$Syxminmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$Syxminmax[which(rownames(EstimParams$Syxminmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("Sxyminmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="Sxy")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="Sxystart")):(which(colnames(EstimParams$mGrid)=="Sxyend"))}
	for (i in colstr){
	    minVal<-EstimParams$Sxyminmax[which(rownames(EstimParams$Sxyminmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$Sxyminmax[which(rownames(EstimParams$Sxyminmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    if (is.element("Sxxminmax",names(EstimParams))){
	colstr<-which(colnames(EstimParams$mGrid)=="Sxx")
	if (length(colstr)==0){colstr<-(which(colnames(EstimParams$mGrid)=="Sxxstart")):(which(colnames(EstimParams$mGrid)=="Sxxend"))}
	for (i in colstr){
	    minVal<-EstimParams$Sxxminmax[which(rownames(EstimParams$Sxxminmax)==colnames(EstimParams$mGrid)[i]),1]
	    maxVal<-EstimParams$Sxxminmax[which(rownames(EstimParams$Sxxminmax)==colnames(EstimParams$mGrid)[i]),2]
	    EstimParams$mGrid[,i]<-(maxVal-minVal)*atan(EstimParams$mGrid[,i])/pi+(minVal+maxVal)/2
	}
    }
    EstimParams$mGrid
}
