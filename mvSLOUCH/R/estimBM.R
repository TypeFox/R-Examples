.bm.estim<-function(dfX,phylTree,Merror=NULL){
## uses the brown function in ouch and also assumes ouch format of phylogenetic tree
    vNAX<-which(is.na(c(t(as.matrix(dfX)[phylTree@term,]))))
    LogLik<- -Inf
    RSS <- Inf
    vX0<-rep(0,ncol(dfX))
    StS<-diag(1,ncol(dfX),ncol(dfX))
    Sxx<-StS

    if (length(vNAX)==0){
	tryCatch({
	    bmEstim<-ouch::brown(data=dfX,tree=phylTree)
	    LogLik<-bmEstim@loglik
	    vX0<-matrix(sapply(bmEstim@theta,function(x){x},simplify=TRUE),ncol=1,nrow=ncol(dfX))
	    StS<-.sym.par(bmEstim@sigma,ncol(dfX))
	    Sxx<-t(chol(StS))
	    V1<-pseudoinverse(phylTree@branch.times%x%StS)
	    data<-c(t(dfX[(phylTree@nnodes-phylTree@nterm+1):phylTree@nnodes,]))
	    RSS<-((data-rep(vX0[,1],phylTree@nterm))%*%V1%*%(data-rep(vX0[,1],phylTree@nterm)))[1,1]
	},error=function(e){print(paste("Error in BM optim ",e))})
    }
    
    StS[which(abs(StS)<1e-15)]<-0
    Sxx[which(abs(Sxx)<1e-15)]<-0
    vX0[which(abs(vX0)<1e-15)]<-0
    vX0<-matrix(vX0,ncol=1,nrow=ncol(dfX))
    if (!(is.null(Merror)||is.na(Merror)||(sum(sum(abs(Merror)))<1e-15))||(length(vNAX)!=0)||(is.nan(LogLik)||is.na(LogLik)||(LogLik< -1000000)||is.infinite(LogLik))){## there is measurement error or something went wrong
	D<-matrix(1,ncol=1,nrow=length(phylTree@term))%x%diag(1,ncol(dfX),ncol(dfX))
	data<-c(t(as.matrix(dfX)[phylTree@term,]))
	if (length(vNAX)>0){data<-data[-vNAX];if(ncol(D)>1){D<-D[-vNAX,]}else{D<-matrix(D[-vNAX,],ncol=1)}}
	optPar<-NA
	##one-diml optimization by Nelder-Mead is unreliable: use optimize however after testing optim performs much better problem with setting boundries for optimize
	tryCatch({
	    optSxx<-optim(par=.sym.unpar(StS),fn=function(parStS,T,Merror,data,D,vNAX){
		LogLik<- 1000000
		StS<-.sym.par(parStS)
		V<-T%x%StS
		if (!is.null(Merror)){V<-V+Merror}
		if (length(vNAX)>0){V<-V[-vNAX,-vNAX]}
		V1<-pseudoinverse(V)
		vX0<-pseudoinverse(t(D)%*%V1%*%D)%*%t(D)%*%V1%*%data
		vX0[which(abs(vX0)<1e-15)]<-0
		Edata<-rep(vX0,nrow(T))		
		vX0<-matrix(vX0,ncol=1,nrow=ncol(dfX))
		if (length(vNAX)>0){Edata<-Edata[-vNAX]}
		tryCatch({
		    LogLik<-(-1)*dmvnorm(data,mean=Edata,sigma=V,log=TRUE)
		},error=function(e){})
		if (is.nan(LogLik)||is.na(LogLik)||(LogLik> 1000000)||is.infinite(LogLik)){LogLik<- 1000000}
		LogLik
	    },T=phylTree@branch.times,Merror=Merror,data=data,D=D,vNAX=vNAX);
	    optPar<-optSxx$par
	    LogLik<- (-1)*optSxx$value
	    },error=function(e){print(paste("Error in BM optim ",e))}	
	)
	if (!is.na(optPar[1])){
	    StS<-.sym.par(optPar,ncol(dfX))
	    Sxx<-t(chol(StS))
	    StS[which(abs(StS)<1e-15)]<-0
    	    Sxx[which(abs(Sxx)<1e-15)]<-0
    	    V<-phylTree@branch.times%x%StS
	    if (!is.null(Merror)){V<-V+Merror}
	    if (length(vNAX)>0){V<-V[-vNAX,-vNAX]}
	    V1<-pseudoinverse(V)

	    vX0<-pseudoinverse(t(D)%*%V1%*%D)%*%t(D)%*%V1%*%data
	    vX0[which(abs(vX0)<1e-15)]<-0
	    vX0<-matrix(vX0,ncol=1,nrow=ncol(dfX))
	    vMean<-rep(vX0,phylTree@nterm)
	    if (length(vNAX)>0){vMean<-vMean[-vNAX]}
	    RSS<- ((data-vMean)%*%V1%*%(data-vMean))[1,1]
	}
    }
    list(vX0=vX0,StS=StS,Sxx=Sxx,LogLik=LogLik,RSS=RSS)
}
                                                                                                                                                                                                                                                                                                                                                                                                                                              