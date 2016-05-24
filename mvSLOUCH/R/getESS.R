.calcESSsim<-function(N=100){
## we will simulate using PhyloSDE here
}


.calcESSanalytical<-function(phyltree,proc.params,evolmodel,Merror=NULL,vNAs=NULL,ESS.method="reg"){
## based on ouch tree - will need to move to ape
    n<-phyltree@nterm
    if (is.null(Merror) || is.na(Merror)) {Merror<-0}
    if (evolmodel=="slouch"){evolmodel<-"mvslouch"}
    if ((evolmodel=="bm") || (evolmodel=="ouch") || (evolmodel=="mvslouch")){
    	res<-.getMVSLphylCovMatrix(phyltree,proc.params,evolmodel,Merror)
	procdim<-res$procdim
	V<-res$mCovPhyl.Merror
    }
    orgV<-V
    if ((!is.null(vNAs)) && (length(vNAs)>0)){V<-V[-vNAs,-vNAs]}
    samplesize<-n*procdim
    if ((!is.null(vNAs)) && (length(vNAs)>0)){samplesize<-samplesize-length(vNAs)}
    full.NA.species<-which(sapply(1:n,function(i,d,vNAs){res<-FALSE;if(length(intersect((((i-1)*d+1):(i*d)),vNAs))==d){res<-TRUE};res},d=procdim,vNAs=vNAs,simplify=TRUE))
    n.noNAs<-n - length(full.NA.species)

## samplesize should equal nV

    logdetV<-log(det(V)) ##determinant(V,logarithm=TRUE)
    sumdiagV<-sum(diag(V))
    sumlogdiagV<-sum(log(diag(V)))
    sumoffdiagV<-sum(V[upper.tri(V,diag=FALSE)])
    sumoffdiagabsV<-sum(abs(V[upper.tri(V,diag=FALSE)]))
    sumabsV<-sum(abs(V))
    sumV<-sum(V)
    nV<-nrow(V)
    corrV<-cov2cor(V)    
    sumoffdiagabscorrV<-sum(abs(corrV[upper.tri(corrV,diag=FALSE)]))
    rhoabsV<-mean(abs(corrV[upper.tri(corrV,diag=FALSE)]))

    ESS.factor<-0
    ESS<-1
    ESS.factor.model.selection<-0
    ESS.model.selection<-1

    ESS.factor.MI<-0
    ESS.MI<-1
    ESS.factor.MI<-1/log(exp(1)+(sumlogdiagV-logdetV)/2)
    ESS.MI<-(samplesize-1)*ESS.factor.MI+1
    if (ESS.factor.MI<0){
	print("Warning : negative MI ESS factor")
	ESS.factor.MI<-0
	ESS.MI<-1
    }
    if (ESS.method=="MI"){ESS.factor<-ESS.factor.MI;ESS<-ESS.MI}

    ESS.factor.mean<-0
    ESS.mean<-1
    tryCatch({
	tmpcorrV<-corrV
	numtry<-0
	while(numtry<3){
	    if (numtry>0){tmpcorrV<-jitter(tmpcorrV,amount=1e-14);print("Jittering for mean ESS")}
    	    numtry<-numtry+1
    	    tryCatch({invR<-pseudoinverse(tmpcorrV);numtry<-5},error=function(e){})
	}	
	if (numtry==3){pseudoinverse(tmpcorrV)} ## generate error
	ESS.mean<-sum(invR)
	ESS.factor.mean<-(ESS.mean-1)/(samplesize-1)
	if (ESS.factor.mean<0){
	    print("Warning : negative mean ESS factor")
	    ESS.factor.mean<-0
	    ESS.mean<-1
	}
    },error=function(e){print("Problem with inverting between species correlation matrix");print(e)})
    if (ESS.method=="mean"){ESS.factor<-ESS.factor.mean;ESS<-ESS.mean;ESS.factor.model.selection<-ESS.factor.mean;ESS.model.selection<-ESS.mean}

    ESS.factor.reg<-0
    ESS.reg<-1
    tryCatch({
    	ESS.reg<-sum(sapply(1:nV,function(i,V){
	    numtry<-0
	    while(numtry<3){
		    if (numtry>0){V[-i,-i]<-jitter(V[-i,-i],amount=1e-14);print("Jittering for reg ESS")}
    		    numtry<-numtry+1
    		    tryCatch({
    			res<-1-V[i,-i]%*%pseudoinverse(V[-i,-i])%*%V[-i,i]/V[i,i]
    			if ((res >=0) || (numtry==3)){numtry<-5}
    		    },error=function(e){})    		
    	    }
    	    if (numtry==3){pseudoinverse(V[-i,-i])} ## generate error
    	    if (res<0){print("Error in pseudoinverse --- negative variance, taking absolute value");print(c(i,res));res<-abs(res)}
    	    res
    	},V=V,simplify=TRUE))
	ESS.factor.reg<-(ESS.reg-1)/(samplesize-1)
	if (ESS.factor.reg<0){
	    print("Warning : negative regression ESS factor")
	    ESS.factor.reg<-0
	    ESS.reg<-1	    
	}
    },error=function(e){print("Problem with calculating multivariare regression ESS.");print(e)})
    if (ESS.method=="reg"){ESS.factor<-ESS.factor.reg;ESS<-ESS.reg;ESS.factor.model.selection<-ESS.factor.reg;ESS.model.selection<-ESS.reg;}


    ESS.factor.mvMI<-0
    ESS.mvMI<-1
    if (procdim>1){
	sumlogdetVj<-0
	ni<-n
	for (i in setdiff(1:n,full.NA.species)){
	    Vj<-orgV[((((i-1)*procdim+1)):(i*procdim)),((((i-1)*procdim+1)):(i*procdim))]
	    logdetVj<-0
	    if ((!is.null(vNAs)) && (length(vNAs)>0)){
		vVNAj<-intersect(vNAs,((((i-1)*procdim+1)):(i*procdim)))
		if ((!is.null(vVNAj)) && (length(vVNAj)>0)){
		    if (length(vVNAj)<procdim){
			vVNAj<-(vVNAj-1)%%procdim+1
			Vj<-Vj[-vVNAj,-vVNAj]
			logdetVj<-log(det(Vj))
		    }else{ni<-ni-1}
		}else{logdetVj<-log(det(Vj))}
	    }else{logdetVj<-log(det(Vj))}
	    sumlogdetVj<-sumlogdetVj+logdetVj
	}
	ESS.mvMI<-1
	ESS.factor.mvMI<-1/log(exp(1)+(sumlogdetVj-logdetV)/2)
	ESS.mvMI<-(n.noNAs-1)*ESS.factor.mvMI+1
        if (ESS.method=="mvMI"){ESS.factor<-ESS.factor.mvMI;ESS<-ESS.mvMI;ESS.factor.model.selection<-ESS.factor.MI;ESS.model.selection<-ESS.MI}
    }

    ESS.factor.mvreg<-0
    ESS.mvreg<-1
    if (procdim>1){
	ni<-n
	tryCatch({
    	    ESS.mvreg<-sum(sapply(setdiff(1:n,full.NA.species),function(i,orgV){
		numtry<-0
		while(numtry<3){
		    species.i<-c((procdim*(i-1)+1):(procdim*i))
		    species.noi<- setdiff(1:samplesize,c((procdim*(i-1)+1):(procdim*i)))
		    if (numtry>0){orgV<-jitter(orgV,amount=1e-14);print("Jittering for mvreg ESS")}
    		    numtry<-numtry+1
    		    tryCatch({    			
    			species.i.noNAs<-setdiff(species.i,vNAs)    			
    			species.noi.noNAs<- setdiff(species.noi,vNAs)
    			res<-diag(1,length(species.i.noNAs),length(species.i.noNAs))-pseudoinverse(orgV[species.i.noNAs,species.i.noNAs])%*%orgV[species.i.noNAs,species.noi.noNAs]%*%pseudoinverse(orgV[species.noi.noNAs,species.noi.noNAs])%*%orgV[species.noi.noNAs,species.i.noNAs]
    			if ((res >=0) || (numtry==3)){numtry<-5}
    		    },error=function(e){})    		
    		}
    		if (numtry==3){pseudoinverse(orgV[species.noi.noNAs,species.noi.noNAs])+pseudoinverse(orgV[species.i.noNAs,species.i.noNAs])} ## generate error
    		if (res<0){print("Error in pseudoinverse --- negative variance, taking absolute value");print(c(i,res));res<-abs(res)}
    		det(res) ## in mvreg we have total variance
    	    },orgV=orgV,simplify=TRUE))
	    ESS.factor.mvreg<-(ESS.mvreg-1)/(n.noNAs-1)
	    if (ESS.factor.mvreg<0){
		print("Warning : negative multivariate regression ESS factor")
		ESS.factor.mvreg<-0
		ESS.mvreg<-1	    
	    }
	},error=function(e){print("Problem with calculating multivariare regression ESS.");print(e)})
	ESS.mvreg<-(n.noNAs-1)*ESS.factor.mvreg+1
        if (ESS.method=="mvreg"){ESS.factor<-ESS.factor.mvreg;ESS<-ESS.mvreg;ESS.factor.model.selection<-ESS.factor.reg;ESS.model.selection<-ESS.reg}
    }

    
    rhon<-matrix(0,procdim,procdim)
    if (procdim==1){rhon[1,1]<-mean(c(corrV[upper.tri(corrV,diag=FALSE)]))}
    else{
	for(i in 1:procdim){## this should be done as an sapply
	    for(j in 1:procdim){
		itraits<-which(((1:nrow(corrV))-1)%%procdim+1 == i)
		jtraits<-which(((1:nrow(corrV))-1)%%procdim+1 == j)
		corrVij<-corrV[itraits,jtraits]
		rhon[i,j]<-mean(c(corrVij[upper.tri(corrVij,diag=FALSE)]))
	    }
	}
    }
    list(ESS.factor=ESS.factor,ESS=ESS,ESS.factor.model.selection=ESS.factor,ESS.model.selection=ESS,rhon=rhon,ESS.coeffs=list(ESS.MI=list(ESS.MI.factor=ESS.factor.MI,ESS.MI=ESS.MI,ESS.factor.mvMI=ESS.factor.mvMI,ESS.mvMI=ESS.mvMI),ESS.mean=list(ESS.mean.factor=ESS.factor.mean,ESS.mean=ESS.mean),ESS.reg=list(ESS.reg.factor=ESS.factor.reg,ESS.reg=ESS.reg),ESS.reg=list(ESS.reg.factor=ESS.factor.reg,ESS.reg=ESS.reg,ESS.mvreg.factor=ESS.factor.mvreg,ESS.mvreg=ESS.mvreg)))
}

.getESScriteria<-function(LogLik,dof,ESS,ESS.factor,rhon,RSS){
    lParamsSummary<-list()
    lParamsSummary$LogLik<-LogLik
    lParamsSummary$dof<- dof
    lParamsSummary$m2loglik<- -2*LogLik
    lParamsSummary$aic<- -2*LogLik+2*lParamsSummary$dof
    lParamsSummary$aic.c<- lParamsSummary$aic +2*lParamsSummary$dof*(lParamsSummary$dof+1)/(ESS-lParamsSummary$dof-1)
    lParamsSummary$sic<- lParamsSummary$m2loglik+log(ESS)*lParamsSummary$dof
    lParamsSummary$bic<-lParamsSummary$m2loglik+lParamsSummary$dof*log(ESS)
    lParamsSummary$RSS<-RSS
    lParamsSummary$ESS<-ESS
    lParamsSummary$ESS.factor<-ESS.factor
    lParamsSummary$rhon<-rhon
    lParamsSummary
}


.getMVSLphylCovMatrix<-function(phyltree,proc.params,evolmodel,Merror=NULL){
    n<-phyltree@nterm
    params<-list()
    if (evolmodel=="slouch"){evolmodel<-"mvslouch"}
    params$EvolModel<-evolmodel 
    phyltree@nodelabels[which(phyltree@nodelabels=="")]<-1:length(which(phyltree@nodelabels==""))
    regimes.times<-sapply(phyltree@epochs,rev,simplify=FALSE)
    regimes<-sapply(regimes.times,function(reg){rep("reg.1",length(reg)-1)},simplify=FALSE) ## we only want covariance regimes make no difference
    root.regime<-regimes[1]
    regimes.types.orig<-c("reg.1")
    regimes.types<-1:1
    regimes<-sapply(regimes,function(vregs,regimes.types.orig){sapply(vregs,function(orgreg,regimes.types.orig){which(regimes.types.orig==orgreg)},regimes.types.orig=regimes.types.orig)},regimes.types.orig=regimes.types.orig,simplify=FALSE)
    if (params$EvolModel=="bm"){kY<-0;kX<-nrow(proc.params$vX0)}
    if (params$EvolModel=="slouch"){kY<-1;kX<-ncol(proc.params$B)}
    if (params$EvolModel=="ouch"){kY<-ncol(proc.params$A);kX<-0}
    if (params$EvolModel=="mvslouch"){kY<-ncol(proc.params$A);kX<-ncol(proc.params$B)}
    params$EstimParams<-list()
    if (!is.element("TerminalLabels",params$EstimParams)){params$EstimParams$TerminalLabels<-phyltree@nodelabels[phyltree@term]}
    lPrecalculates<-mvSLOUCH::.calculate.Tree.dists(phyltree,UserTermLabels=params$EstimParams$TerminalLabels)
    modelParams<-vector("list",4); names(modelParams)<-c("regimeTimes","regimes","regimeTypes","Merror")
    modelParams$regimeTimes<-regimes.times; modelParams$regimes<-regimes ; modelParams$regimeTypes<-regimes.types
    if (is.null(Merror) || is.na(Merror)) {Merror<-0}
    modelParams$Merror<-.createCovariancematrix(Merror,n,kX+kY,NULL,"measurement error") ## Here we will need to change package to GLSME
    modelParams<-c(modelParams,proc.params)
    bCalcA<-TRUE ; lexptcalc<-TRUE; if (evolmodel=="bm"){bCalcA<-FALSE;lexptcalc<-FALSE;}
    modelParams$precalcMatrices<-mvSLOUCH::.decompEigenA.S(modelParams,lPrecalculates,NULL,list(bCalcA=bCalcA,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=lexptcalc,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    V<-mvSLOUCH::.calc.phyl.cov(lPrecalculates$mTreeDist,lPrecalculates$mSpecDist[nrow(lPrecalculates$mSpecDist),],lPrecalculates$mAncestorTimes,lPrecalculates$vSpeciesPairs,evolmodel,modelParams)
    list(mCovPhyl=V,mCovPhyl.Merror=V+modelParams$Merror,procdim=kY+kX)
}

