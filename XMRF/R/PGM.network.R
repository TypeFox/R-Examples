PGM.network <-
function(X,R=max(X),nlams=10,lambda=NULL,parallel=TRUE,ncores=ncores){

	ghat = c()
	if(!is.null(lambda)){
		ghat = array(0,dim=c(nrow(X),nrow(X),length(lambda)))
		nlams = length(lambda)
	}
	else{
		lmax = lambdaMax(t(X))
		lambda = exp(seq(log(lmax),log(0.01*lmax),l=nlams))
		ghat = array(0,dim=c(nrow(X),nrow(X),length(lambda)))
	}
	
	#lmax = lambdaMax(t(X))
	#lams = exp(seq(log(lmax),log(.0001),l=nlams));
	
	  #-if(nlams>0){
	  #		ghat = array(0,dim=c(nrow(X),nrow(X),length(lams)))
	  # }
	wrapper <- function(i){
		
		#print(i)	
		#fit = PGM.path.neighborhood(t(X[-i,]),X[i,],nlams,lmax,0)
		fit = PGM.path.neighborhood(t(X[-i,]),X[i,],nlams=nlams,startb=0,lambda=lambda)
		fit$beta=as.matrix(fit$Bmat)
		if(i==1){
			ghat[i,2:nrow(X),]=fit$beta
		}else if(i==nrow(X)){
			ghat[i,1:(nrow(X)-1),]=fit$beta

		}else{
			ghat[i,1:(i-1),]=fit$beta[1:(i-1),]
			ghat[i,(i+1):nrow(X),]=fit$beta[i:nrow(fit$beta),]	
		}
		return(ghat[i,,])
	}
	
	ghat2=c()
	if(parallel){
		###library(multicore)
		ghat2=parallel::mclapply(1:nrow(X),wrapper,mc.cores=ncores)	
	}else{
		ghat2=lapply(1:nrow(X),wrapper)	
	}
	for(i in 1:nrow(X)){
		ghat[i,,]=ghat2[[i]]
	}
	
	ghat=lapply(1:nlams,function(r){return(ghat[,,r])})
	
	return(ghat)
}
