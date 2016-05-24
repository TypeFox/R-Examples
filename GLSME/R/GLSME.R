GLSME<-function(y,D,Vt,Ve,Vd,Vu,EstimateVariance=c(TRUE,TRUE),CenterPredictor=TRUE,InitialGuess=NULL,eps=0.001,MaxIter=50,MaxIterVar=50,epsVar=0.001,OutputType="short",Vttype=NULL,Vetype=NULL,Vdtype=NULL,Vutype=NULL,ED=NULL,EDtype="SingleValue"){
## get dimensions and check whether everything is in correct matrix form
    if (is.vector(y)){y<-matrix(y,nrow=length(y),ncol=1)}
    n<-nrow(y)
    if (is.vector(D)){D<-matrix(D,nrow=n,ncol=1)}
    m<-ncol(D)    
## ---------------------------------------------------------------------------
    
## calculate covariance matrices from information provided by the user
    lVt<-.createCovariancematrix(Vt,n,1,Vttype,"Vt")
    Vt<-lVt[[1]];Vttype<-lVt[[2]]
    lVe<-.createCovariancematrix(Ve,n,1,Vetype,"Ve")
    Ve<-lVe[[1]];Vetype<-lVe[[2]]
    lVd<-.createCovariancematrix(Vd,n,m,Vdtype,"Vd")
    Vd<-lVd[[1]];Vdtype<-lVd[[2]]
    lVu<-.createCovariancematrix(Vu,n,m,Vutype,"Vu")
    Vu<-lVu[[1]];Vutype<-lVu[[2]]

    Vt[which(abs(Vt)<1e-13)]<-0
    Vd[which(abs(Vd)<1e-13)]<-0
    Ve[which(abs(Ve)<1e-13)]<-0
    Vu[which(abs(Vu)<1e-13)]<-0
## ---------------------------------------------------------------------------

    predVarKnown<-TRUE
    if (length(EstimateVariance[EstimateVariance[2:length(EstimateVariance)]])){predVarKnown<-FALSE}    
    if ((length(EstimateVariance)==2)&&(m>1)){EstimateVariance<-c(EstimateVariance[1],rep(EstimateVariance[2],m))}
    fixedrow<-c() ## find fixed effects, i.e. 0 rows in covariance matrix Vd
    for (i in 1:nrow(Vd)){if (base::sum(abs(Vd[i,]))<1e-13){fixedrow<-c(fixedrow,i)}}
## calculate expectation of design and variance components of predictors
    if (!is.null(ED)){ ## the user provided some information about the mean of the design
	mED<-matrix(NA,ncol=m,nrow=n)
	if (((EDtype=="Matrix")||(is.null(EDtype)))&&(((is.matrix(ED))&&(ncol(ED)==m)&&(nrow(ED)==n))||((m==1)&&(length(c(ED))==n)))){mED<-ED} ## user gave the full expectation
	else{## user provided some structured expectation
	    if (EDtype=="SingleValue"){mED<-matrix(ED,nrow=n,ncol=m)} ## all design has the same mean
	    if (EDtype=="Vector"){## each predictor has its own mean
		mED<-matrix(ED,nrow=n,ncol=m,byrow=TRUE)
		for (i in 1:m){if (ED[i]=="F"){mED[,i]<-D[,i]}}## if fixed effect then 0 variance mean equals observed
		mED<-as.numeric(mED) ## correct to be numeric if "F"s were left trailing
	    }
	}
	ED<-mED
	if (!predVarKnown){## we know Vd up to a constant(s)
	    PredVarEstimates<-c()
	    if ((Vdtype=="Vector")||(Vdtype=="MatrixList")){## if we can treat all predictor variables separately
		PredVarEstimates<-rep(NA,m) 
		for (i in 1:m){ ##  for each predictor
		    if (base::sum(abs(Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]))>1e-13){ ## check if this is not a fixed effect
			if (EstimateVariance[i+1]){
			    options(warn=-1)
			    PredVarEstimates[i]<-exp(optim(par=0,fn=function(vx,Vd,Vu,x,EDx){ ## maximum likelihood estimation
				vx<-exp(vx) ## logged so estimate will be forced to be positive
				(-1)*dmvnorm(x,mean=EDx,sigma=vx*Vd+Vu,log=TRUE)
			    },Vd=Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)],Vu=Vu[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)],x=D[,i],EDx=ED[,i])$par)		    
			    options(warn=1)
			    Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]<-PredVarEstimates[i]*Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)] ## update Vd
			}else{PredVarEstimates[i]<-1}
		    }else{PredVarEstimates[i]<-0;Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]<-0}
		}
	    }else{ ## we cannot treat the predictors separately
	        if (base::sum(abs(Vd))>1e-13){ ## we check if by chance not all of them are fixed effects
		    options(warn=-1)
		    if (length(fixedrow)>0){Vduse<-Vd[-fixedrow,-fixedrow];Vuuse=Vu[-fixedrow,-fixedrow];xuse=c(D)[-fixedrow];EDxuse=c(ED)[-fixedrow]}
		    else{Vduse<-Vd;Vuuse=Vu;xuse=c(D);EDxuse=c(ED)}
		    PredVarEstimates<-exp(optim(par=0,fn=function(vx,Vd,Vu,x,EDx){
			vx<-exp(vx)
			(-1)*dmvnorm(x,mean=EDx,sigma=vx*Vd+Vu,log=TRUE)
		    },Vd=Vduse,Vu=Vuuse,x=xuse,EDx=EDxuse)$par)		    
		    options(warn=1)
		    Vd<-PredVarEstimates*Vd ## update Vd
		}else{PredVarEstimates<-0;Vd[1:nrow(Vd),1:ncol(Vd)]<-0}
	    }	    
	}
    }
    else{ ## we need to estimate ED
	tmpVd<-Vd
	if (!predVarKnown){
    	    if ((Vdtype=="Vector")||(Vdtype=="MatrixList")){varconsts<-rep(1,m)}else{varconsts<-1}
    	    tmpvarconsts<-varconsts*100+1000
    	    iter<-1
	    while((iter<MaxIterVar)&&((base::sum((varconsts-tmpvarconsts)^2))>epsVar)){ ## iterative search to find ED and variance constants
		tmpvarconsts<-varconsts
		tmpD<-c(D)
		tmpED<-rep(NA,m*n)
		for (i in 1:m){## estimate the mean separately for each variable
		    tmpVdu1i<-pseudoinverse(tmpVd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]+Vu[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)],0.000001)
    		    tmpED[((i-1)*n+1):(i*n)]<-tmpD[((i-1)*n+1):(i*n)]%*%tmpVdu1i%*%rep(1,n)/(base::sum(tmpVdu1i)) ## current best guess of ED[,i]
		}		    		

    		if(length(fixedrow)>0){tmpED[fixedrow]<-tmpD[fixedrow]} ## correct for known fixed effects
    		tmpED<-matrix(tmpED,ncol=m,nrow=n,byrow=FALSE)    
    		## and now conditional on ED estimate the variance constants 
    	        if ((Vdtype=="SingleVector")||(Vdtype=="MatrixList")){ ## if we can treat all predictor variables separately
		    for (i in 1:m){##  for each predictor
			if (base::sum(abs(Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]))>1e-13){## check if not fixed effect
			    if (EstimateVariance[i+1]){
				options(warn=-1)
				varconsts[i]<-exp(optim(par=0,fn=function(vx,Vd,Vu,x,EDx){ ## maximum likelihood estimation
				    vx<-exp(vx) ## logged to force to be positive
				    (-1)*dmvnorm(x,mean=EDx,sigma=vx*Vd+Vu,log=TRUE)
				},Vd=Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)],Vu=Vu[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)],x=D[,i],EDx=tmpED[,i])$par)		    
				options(warn=1)
				tmpVd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]<-varconsts[i]*Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]
			    }else{varconsts[i]<-1}
			}else{varconsts[i]<-0;Vd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]<-0;tmpVd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]<-0}			
		    }
		}else{ ## we cannot treat the predictors separately
		    options(warn=-1)
		    if (length(fixedrow)>0){Vduse<-Vd[-fixedrow,-fixedrow];Vuuse=Vu[-fixedrow,-fixedrow];xuse=c(D)[-fixedrow];EDxuse=c(tmpED)[-fixedrow]}
		    else{Vduse<-Vd;Vuuse=Vu;xuse=c(D);EDxuse=c(tmpED)}
		    varconsts<-exp(optim(par=0,fn=function(vx,Vd,Vu,x,EDx){
			vx<-exp(vx) ## logged to force to be positive
			(-1)*dmvnorm(x,mean=EDx,sigma=vx*Vd+Vu,log=TRUE)
		    },Vd=Vduse,Vu=Vuuse,x=xuse,EDx=EDxuse)$par)		    
		    options(warn=1)
		    tmpVd<-varconsts*Vd
		}
		iter<-iter+1
    	    }	
	    orgVd<-Vd
	    PredVarEstimates<-varconsts	
    	    Vd<-tmpVd ## update Vd
	}
	tmpD<-c(D)
	Vdu1<-pseudoinverse(Vd+Vu,0.000001)
    	mn<-m*n
	tmpED<-rep(NA,m*n)
	for (i in 1:m){## estimate the mean separately for each variable
	    tmpVdu1i<-pseudoinverse(tmpVd[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]+Vu[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)],0.000001)
	    tmpED[((i-1)*n+1):(i*n)]<-tmpD[((i-1)*n+1):(i*n)]%*%tmpVdu1i%*%rep(1,n)/(base::sum(tmpVdu1i)) ## current best guess of ED[,i]
	}		    			
	if (length(fixedrow)>0){tmpED[fixedrow]<-0} ## correct for known fixed effects
    	tmpED<-matrix(tmpED,ncol=m,nrow=n,byrow=FALSE)            
	ED<-tmpED
    }     
    orgD<-D
    if (CenterPredictor){
	D<-D-ED ## center the variables
	if (length(fixedrow)>0){
	    ED[fixedrow]<-D[fixedrow]
    	    ED[-fixedrow]<-0
        }
    }
## ---------------------------------------------------------------------------    
    ## do iterative GLS estimation for regression coefficients
    ## setup initial starting point
    if (is.null(InitialGuess)){InitialGuess<-pseudoinverse(t(D)%*%D,0.000001)%*%t(D)%*%y}
    if ((EstimateVariance[1])&&(length(InitialGuess)==m)){InitialGuess<-c(InitialGuess,1)}
    GLSestim<-matrix(InitialGuess,ncol=1)
    tmpGLSestim<-InitialGuess*100+1000
    tmpVt<-Vt
    Vud1<-pseudoinverse(Vu+Vd,0.000001)

    ## calculate Var(Ub | D)
    VarUbcD<-Reduce('+',sapply(1:m^2,function(rk,m,n,Vucd,b){
	    r<-(rk-1)%/%m+1
	    k<-(rk-1)%%m+1
	    b[r]*b[k]*Vucd[((r-1)*n+1):(r*n),((k-1)*n+1):(k*n)]
    },m=m,n=n,Vucd=Vu-Vu%*%Vud1%*%Vu,b=GLSestim,simplify=F),init=matrix(0,n,n))

    iter<-1
    while((iter<MaxIter)&&((base::sum((GLSestim-tmpGLSestim)^2))>eps)){ ## main search loop
	tmpGLSestim<-GLSestim
	if (EstimateVariance[1]){vy0<-GLSestim[m+1,1];tmpVt<-vy0*Vt} ## update current guess of Vt
	tmpVte<-tmpVt+Ve
	
	## normal GLS calculation
    	V<-tmpVte+VarUbcD
        V1<-pseudoinverse(V,0.000001)
        tDV1<-t(D)%*%V1
	invtDV1DtDV1<-pseudoinverse(tDV1%*%D,0.0000001)%*%tDV1
        GLSestim<-invtDV1DtDV1%*%y
	
	## calculate Var(Ub | D)
	VarUbcD<-Reduce('+',sapply(1:m^2,function(rk,m,n,Vucd,b){
	    r<-(rk-1)%%m+1
	    k<-(rk-1)%/%m+1
	    b[r]*b[k]*Vucd[((r-1)*n+1):(r*n),((k-1)*n+1):(k*n)]
	},m=m,n=n,Vucd=Vu-Vu%*%Vud1%*%Vu,b=GLSestim,simplify=F),init=matrix(0,n,n))

	if (EstimateVariance[1]){ ## get estimate of variance constant for response	    

	    ## calculate Cov(y,D)
	    mCovYD<-Reduce('+',sapply(1:m,function(r,m,n,Vd,b){
		b[r]*Vd[((r-1)*n+1):(r*n),]
	    },m=m,n=n,Vd=Vd,b=GLSestim,simplify=F),init=matrix(0,n,m*n))
    	    condPartEy<-mCovYD%*%Vud1%*%(c(D-ED))  ## E[Y|Do]

	    Vpart<-Ve+VarUbcD
	    vy0<-log(vy0) ## logged to force search to be positive 	    
	    
	    options(warn=-1)
	    GLSestim<-rbind(GLSestim,exp(optim(par=vy0,fn=function(vy,Vpart,Vt,y,D,b,condPartEy,Vud){ ## maximum likelihood estimation
		vy<-exp(vy)  ## logged to force search to be positive 	  
	        ## for this to work we need to know y's expectation | Do but because we have D observed with error
                ## we do not know it actually so  we guess it by making a random draw of the error matrix
                ## of the design, this actually seems to work if the error is not too large
            	DeGuess<-matrix(c(rmvnorm(1,mean=rep(0,nrow(Vud)),sigma=Vud)),nrow=n,byrow=FALSE)
    	        muyD<- (D-DeGuess)%*%b+condPartEy
	        (-1)*dmvnorm(c(y),mean=muyD,sigma=vy*Vt+Vpart,log=TRUE)		  
	    },Vpart=Vpart,Vt=Vt,y=y,D=D,b=GLSestim[1:m,1],condPartEy=condPartEy,Vud=Vu+Vd)$par))
	    options(warn=1)
	}
	print(paste("Finished running iteration ",iter," of iterated GLS. Current estimate : ",sep=""))
	if (!is.null(colnames(D))){ ## if the variables had names, then name appropriate rows and columns of the output
	    if (EstimateVariance[1]){rownames(GLSestim)<-c(colnames(D),"ResponseVariance")}
	    else{rownames(GLSestim)<-colnames(D)}
	}
	print(GLSestim)
	iter<-iter+1
    }
## ---------------------------------------------------------------------------    

    ## Calculate final parameter estimates

    if (EstimateVariance[1]){RespVarEst<-GLSestim[m+1,1];Vt<-RespVarEst*Vt;names(RespVarEst)<-NULL}
    Vte<-Vt+Ve

    V<-Vte+VarUbcD
    V1<-pseudoinverse(V,0.000001)
    tDV1<-t(D)%*%V1
    invtDV1D<-pseudoinverse(tDV1%*%D,0.0000001)
    invtDV1DtDV1<-invtDV1D%*%tDV1
    GLSestim<-invtDV1DtDV1%*%y
    EUcD<-matrix(Vu%*%pseudoinverse(Vu+Vd,0.000001)%*%(c(D)-c(ED)),ncol=ncol(D),nrow=nrow(D),byrow=FALSE)
    K<-diag(1,m,m)-invtDV1DtDV1%*%EUcD
    K1<-pseudoinverse(K,0.000001)
    BiasCorrGLS<-K1%*%GLSestim[1:m,1]
## ---------------------------------------------------------------------------

## calculate estimates' covariance end standard error estimates
    covarerrorbiasGLS<-t(K1)%*%invtDV1D%*%K1
    errorGLS<-matrix(sqrt(diag(invtDV1D)),ncol=1)
    errorBiasGLS<-matrix(sqrt(diag(covarerrorbiasGLS)),ncol=1)
## ---------------------------------------------------------------------------    

## calculate R2
    Ey<-(t(y)%*%V1%*%rep(1,n)/(base::sum(V1)))[1,1]
    R2<-1-(t(y-(D-EUcD)%*%GLSestim)%*%V1%*%(y-(D-EUcD)%*%GLSestim))/(t(y-Ey)%*%V1%*%(y-Ey))
    biasR2<-1-(t(y-(D-EUcD)%*%BiasCorrGLS)%*%V1%*%(y-(D-EUcD)%*%BiasCorrGLS))/(t(y-Ey)%*%V1%*%(y-Ey))
## calculate loglikelihood
    LogLik<- dmvnorm(c(y),mean=c((D-EUcD)%*%GLSestim),sigma=V,log=TRUE)   
    LogLikBias<- dmvnorm(c(y),mean=c((D-EUcD)%*%BiasCorrGLS),sigma=V,log=TRUE) 
## ---------------------------------------------------------------------------    
    if (!is.null(colnames(D))){ ## if the variables had names, then name appropriate rows and columns of the output
        rownames(GLSestim)<-colnames(D)
	rownames(BiasCorrGLS)<-colnames(D)
	rownames(errorGLS)<-colnames(D)
	rownames(errorBiasGLS)<-colnames(D)
	rownames(K)<-colnames(D)
	colnames(K)<-colnames(D)
	rownames(invtDV1D)<-colnames(D)
	colnames(invtDV1D)<-colnames(D)		
	rownames(covarerrorbiasGLS)<-colnames(D)
	colnames(covarerrorbiasGLS)<-colnames(D)
	if((!predVarKnown)&&((Vdtype=="SingleVector")||(Vdtype=="MatrixList"))){names(PredVarEstimates)<-colnames(D)}
    }
## ---------------------------------------------------------------------------        
    ## create output list
    output<-list(GLSestimate=GLSestim,errorGLSestim=errorGLS,BiasCorrectedGLSestimate=BiasCorrGLS,errorBiasCorrectedGLSestim=errorBiasGLS,K=K,R2=R2,R2BiasCorrectedModel=biasR2,LogLikelihood=LogLik,LogLikelihoodBiasCorrectedModel=LogLikBias)
    if (!predVarKnown){output<-c(output,list("PredictorVarianceConstantEstimate"=PredVarEstimates))}
    if (EstimateVariance[1]){output<-c(output,list("ResponseVarianceConstantEstimate"=RespVarEst))}
    if (OutputType=="long"){
	output<-c(output,list("CovarianceGLSestimate"=invtDV1D,"CovarianceBiasCorrectedGLSestimate"=covarerrorbiasGLS))
	output<-c(output,list("response"=y,"design"=D,"Vt"=Vt,"Ve"=Ve,"Vd"=Vd,"Vu"=Vu))
    }
    output
}

.createCovariancematrix<-function(covmat,n,m,matrixtype,whichcov){
## function creates the covariance matrix from the data provided by the user according to the provided matrix type
    mCov<-matrix(NA,ncol=n*m,nrow=n*m)
    if (is.null(matrixtype)){
	matrixOK<-FALSE
	if ((is.vector(covmat))&&(length(covmat)==1)){matrixtype<-"SingleValue";matrixOK<-TRUE}
	if ((is.vector(covmat))&&(length(covmat)==m)&&(m>1)){matrixtype<-"Vector";matrixOK<-TRUE;}
	if ((is.list(covmat))&&(length(covmat)==m)&&(m>1)){matrixtype<-"MatrixList";matrixOK<-TRUE}
	if ((is.matrix(covmat))&&(nrow(covmat)==n*m)&&(ncol(covmat)==n*m)&&(n*m>1)){matrixtype<-"Matrix";matrixOK<-TRUE}
	if (!matrixOK){print(paste("Could not setup ",whichcov," matrix. Please check input or provide matrix type.",sep=""));stop()}
    }
    
    if ((matrixtype=="Matrix")&&(is.matrix(covmat))&&(ncol(covmat)==n*m)&&(nrow(covmat)==n*m)){mCov<-covmat}
    else{
        if (matrixtype=="VectorList"){matrixtype<-"MatrixList"}
	if (matrixtype=="SingleValue"){mCov[1:(n*m),1:(n*m)]<-0;diag(mCov)<-covmat}
	if (matrixtype=="Vector"){mCov[1:(n*m),1:(n*m)]<-0;covmat[which(covmat=="F")]<-0;covmat<-as.numeric(covmat);diagvect<-c();for (i in 1:m){diagvect<-c(diagvect,rep(covmat[i],n))};diag(mCov)<-diagvect}
	if (matrixtype=="CorrelatedPredictors"){mCov<-covmat%x%diag(1,n,n)}
	if (matrixtype=="MatrixList"){
	    mCov[1:(n*m),1:(n*m)]<-0;
	    for (i in 1:m){
		if ((is.character(covmat[[i]]))&&(covmat[[i]]=="F")){covmat[[i]]<-0};
		if (is.vector(covmat[[i]])){
		    if((length(covmat[[i]])==n)||(length(covmat[[i]])==1)){
			vcov<-covmat[[i]]
			covmat[[i]]<-matrix(0,n,n)
			diag(covmat[[i]])<-vcov
		    }else{print(paste("Could not setup covariance of predictor ",i," of ",whichcov," matrix. Please check input or provide matrix type.",sep=""));stop()}		
		}
		mCov[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)]<-covmat[[i]]
	    }
	}
	if (matrixtype=="BM"){if (!is.matrix(covmat[[1]])){covmat[[1]]<-diag(covmat[[1]],m,m)};mCov<-covmat[[1]]%x%covmat[[2]]}	    
    }    
    mCov<- (mCov+t(mCov))/2
    EigVals<-eigen(mCov)$values
    if (length(which(EigVals<0)>0)){print(paste("Warning : ",whichcov," matrix has negative eigenvalues!!",sep=""))}
    if (length(which(is.complex(EigVals))>0)){print(paste("Warning : ",whichcov," matrix has complex eigenvalues!!",sep=""))}
    list(mCov,matrixtype)
}





GLSME.predict<-function(xo,glsme.estimate,vy,vx,alpha=0.95){

    yp.biascorr<-xo%*%glsme.estimate$BiasCorrectedGLSestimate
    yp.biasuncorr<-xo%*%glsme.estimate$GLSestimate

    mse.biascorr.estimate<- t(xo)%*%(glsme.estimate$errorBiasCorrectedGLSestim)%*%xo 
			+t(glsme.estimate$BiasCorrectedGLSestimate)%*%vx%*%glsme.estimate$BiasCorrectedGLSestimate+vy

    mse.biasuncorr.estimate<- t(xo)%*%
			(glsme.estimate$CovarianceGLSestimate+
			(matrix((diag(1,length(xo),length(xo))-glsme.estimate$K)%*%glsme.estimate$GLSestimate,ncol=1))%*%
			(matrix((diag(1,length(xo),length(xo))-glsme.estimate$K)%*%glsme.estimate$GLSestimate,nrow=1)))%*%xo 
			+t(glsme.estimate$GLSestimate)%*%vx%*%glsme.estimate$GLSestimate+vy


    n<-length(glsme.estimate$response)
    a<-sqrt(1+1/n)*qt(1-alpha/2,df=n-1)*(vy+t(glsme.estimate$BiasCorrectedGLSestimate)%*%vx%*%glsme.estimate$BiasCorrectedGLSestimate)
    CI.biascorr<-c(yp.biascorr-a,yp.biascorr+a)
    a<-sqrt(1+1/n)*qt(1-alpha/2,df=n-1)*(vy+t(glsme.estimate$GLSestimate)%*%vx%*%glsme.estimate$GLSestimate)
    CI.biasuncorr<-c(yp.biasuncorr-a,yp.biasuncorr+a)
    list(BiasCorr=list(prediction=yp.biascorr,MSE=mse.biascorr.estimate,CI=CI.biascorr,alpha=alpha),BiasUncorr=list(prediction=yp.biasuncorr,MSE=mse.biasuncorr.estimate,CI=CI.biasuncorr,alpha=alpha))
}
