PGM <-
function(X,method="PGM", R=max(X),stability="bootstrap",N=100,beta=0.05,lmin=0.00001,lambda.path=NULL,nlams=20,ncores=4,parallel=F,sth=0.8){
	
	if(is.null(lambda.path) ){
		lmax = lambdaMax(t(X))
 		lambda.path = exp(seq(log(lmax),log(lmin),l=nlams));
	}
	
	if(stability == "bootstrap"){
		b = ncol(X)
		replaceF = TRUE
		
		if(length(lambda.path) > 1){
			cat("Error. Bootstrap allows only one regularization parameter.\n")
			cat("Use STAR stability for the whole regularization.\n")
		}
	}
	else{
		#	b = min(c(10*sqrt(ncol(X)), 0.8*ncol(X))) 
		#	replaceF = FALSE
		b = ncol(X)
		replaceF = TRUE
		
		if(length(lambda.path) == 1){
			cat("Error. STAR allows only for regularization path.\n")
			cat("Use Bootstrap stability for the one regularization parameter.\n")
		}
	}
	
	#	b = min(c(10*sqrt(ncol(X)), 0.8*ncol(X))) 
	ghat=list()
	ghat.path=list()
	ghat.path$path=vector("list",length(lambda.path))
	ghat.path$D = list()
	v=c()
	
	for(i in 1:N){
		cat(paste("PGM: Conducting sampling ... in progress: ", floor(100*(i/N)), "%", collapse=""),"\r")
		flush.console()
		index = sample(1:ncol(X),b,replace=F)
		#tmp=glmpois(X[,index],lambda.path[j],parallel=parallel,warmStart=warmStart,nCpus=nCpus)
		ghat.path$raw = PGM.network(X[,index],R,nlams=length(lambda.path),lambda=lambda.path, parallel=parallel, ncores=ncores)
		
		for(j in 1:length(lambda.path)){
			tmp=ghat.path$raw[[j]]
			tmp[abs(tmp)<1e-06]=0
			tmp[abs(tmp)>1e-06]=1
			diag(tmp)=0
			if(is.null(ghat.path$path[[j]])){
				ghat.path$path[[j]]=tmp;
			}else{
				ghat.path$path[[j]]=ghat.path$path[[j]]+tmp	
			}
		}
	}
		
	for(i in 1:length(lambda.path)){
		if(stability=="bootstrap"){
			D=ghat.path$path[[i]]
			D=D/N
			ghat.path$D[[i]] <- D
			v = NULL
		}
		else{
			D=ghat.path$path[[i]]
			D=D/N
			D=2*D*(1-D)
			ghat.path$D[[i]] <- D
			v=c(v,mean(D[upper.tri(D)]))
		}
	}
					
		
	v=cummax(v)
	ghat$v=v
	ghat$lambda.path = lambda.path
	ghat$D = ghat.path$D
	
	if(stability=="bootstrap"){
		ghat$opt.index = 1
		ghat$opt.lambda = lambda.path[1]
		ghat$network = PGM.network(X,R,nlams,lambda=lambda.path,parallel=parallel,ncores=ncores)
		ghat$network = lapply(1:length(ghat$network),function(i) WeightNet(ghat$network[[i]], ghat$D[[i]], thw=sth))		
	}
	else{
		ghat$opt.index = which(v==max(v[v<beta]))
		ghat$opt.index = max(ghat$opt.index)
		ghat$opt.lambda = lambda.path[which(v==max(v[v<beta]))]
		ghat$network= PGM.network(X,R,nlams,lambda=lambda.path,parallel=parallel,ncores=ncores)
	}
	
	ghat$call <- match.call()
	cat("\nPGM Completed. \n")
	class(ghat) <- "GMS"
	
	return(ghat)
}
