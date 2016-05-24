LPGM <-
function(X, method="LPGM", stability="bootstrap", N=100, beta=0.1, lmin = 0.001, nlams=20, lambda.path=NULL, parallel=TRUE,nCpus=4, sym=TRUE, th=0, sth=0.8){
	
	if(is.null(lambda.path) ){
		lmax = lambdaMax(t(X))
 		lambda.path = exp(seq(log(lmax),log(lmin),l=nlams));
		#lambda.path = exp(seq(log(lmax),log(lmin*lmax),l=nlams))
    ## cat("here\n")
	}
	
	
	if(parallel == T){
		if(stability=="bootstrap"){
			b = ncol(X)
			replaceF = TRUE
			
			if(length(lambda.path) > 1){
				cat("Error. Bootstrap allows only one regularization parameter.\n")
				cat("Use STAR stability for the whole regularization.\n")
				ghat = NULL
				return(ghat)
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
				ghat = NULL
				return(ghat)
			}
		}
		ghat=list()
		ghat.path=list()
		ghat.path$path=vector("list",length(lambda.path))
		ghat.path$D = list()
		v=c()
		
		for(i in 1:N){
			cat(paste(method, "::: Conducting sampling ... in progress: ", floor(100*(i/N)), "%", collapse=""),"\r")
			flush.console()
			
			glmpois.good <- 1
			
			while(glmpois.good){
				# Make sure sample with no gene with all zero values
				good <- 1
				while(good){
					index = sample(1:ncol(X),b,replace=replaceF)
					if(sum(apply(X[,index], 1, function(x) length(unique(x))==1))==0){
						good <- 0
					}
				}
				
				tryCatch(
						{
							ghat.path$raw = LPGM.network(X[,index],lambda=lambda.path, parallel=parallel,nCpus=nCpus,sym=sym,th=th)
							glmpois.good <- 0
						},
						error = function(e) {
							cat("LPGM.network returns empty model. Try again.")
						}
				)
			}
			
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
			ghat$network = LPGM.network(X,lambda=lambda.path,parallel=parallel,nCpus=nCpus,sym=sym,th=th)
			ghat$network = lapply(1:length(ghat$network), function(i) WeightNet(ghat$network[[i]], ghat$D[[i]], thw=sth))
		}
		else{
			ghat$opt.index = which(v==max(v[v<beta]))
			ghat$opt.index = max(ghat$opt.index)
			ghat$opt.lambda = lambda.path[which(v==max(v[v<beta]))]
			ghat$network= LPGM.network(X,lambda=lambda.path,parallel=parallel,nCpus=nCpus,sym=sym,th=th)
		}
			
		ghat$call <- match.call()
		cat(paste("\n", method, " Completed.", "\n", sep=""))
		class(ghat) <- "GMS"
		##return(ghat)
	}
	#
	#
	
	if(parallel == F){
		if(stability=="bootstrap"){
			b = ncol(X)
			replaceF = TRUE
			
			if(length(lambda.path) > 1){
				cat("Error. Bootstrap allows only one regularization parameter.\n")
				cat("Use STAR stability for the whole regularization.\n")
				ghat = NULL
				return(ghat)
			}
		}else{
		#	b = min(c(10*sqrt(ncol(X)), 0.8*ncol(X))) 
			#	replaceF = FALSE
			b = ncol(X)
			replaceF = TRUE
			
			if(length(lambda.path) == 1){
				cat("Error. STAR allows only for regularization path.\n")
				cat("Use Bootstrap stability for the one regularization parameter.\n")
				ghat = NULL
				return(ghat)
			}
		}
		
		ghat=list()
		ghat.path=list()
		ghat.path$D = list()
		v=c()
		
		for( j in 1:length(lambda.path)){
			#	cat(paste(method, ": Conducting sampling ... in progress: ", floor(100*(j/length(lambda.path))), "%", collapse=""),"\r")
			#	flush.console()
			D=matrix(0,nrow=nrow(X),ncol=nrow(X))
			
			for(i in 1:N){

				cat(paste(method, ": Conducting sampling ... in progress: ", floor(100*((i+(j-1)*N)/(N*length(lambda.path)))), "%", collapse=""),"\r")
				flush.console()

				glmpois.good <- 1
				
				while(glmpois.good){
					# Make sure sample with no gene with all zero values
					good <- 1
					while(good){
						index = sample(1:ncol(X),b,replace=replaceF)
						if(sum(apply(X[,index], 1, function(x) length(unique(x))==1))==0){
							good <- 0
						}
					}
					
					tryCatch(
						{
							tmp = LPGM.network(X[,index],lambda=lambda.path[j],parallel=F,sym=sym,th=th)[[1]]
							glmpois.good <- 0
						},
						error = function(e) {
							cat("LPGM.network returns empty model. Try again.\n")
						}
					)
				} 
				
				tmp[abs(tmp)<1e-06]=0
				tmp[abs(tmp)>1e-06]=1
				diag(tmp) = 0
				D=D+tmp
			}
			
			if(stability=="bootstrap"){
				D=D/N
				ghat.path$D[[j]] <- D
				v=c(v,mean(D[upper.tri(D)]))
			}
			else{
				D=D/N
				D=2*D*(1-D)
				ghat.path$D[[j]] <- D
				v=c(v,mean(D[upper.tri(D)]))
			}
			
			#- D=D/N
			#- D=2*D*(1-D)
			#- v=c(v,mean(D[upper.tri(D)]))			
		}
		
		v=cummax(v)
		ghat$v=v
		ghat$lambda.path = lambda.path
		ghat$D = ghat.path$D
		
		if(stability=="bootstrap"){
			ghat$opt.index = 1
			ghat$opt.lambda = lambda.path[1]	
			ghat$network = LPGM.network(X,lambda=lambda.path,parallel=parallel,nCpus=nCpus,sym=sym,th=th)
			ghat$network = lapply(1:length(ghat$network), function(i) WeightNet(ghat$network[[i]], ghat$D[[i]], thw=sth))
		}
		else{
			ghat$opt.index = which(v==max(v[v<beta]))
			ghat$opt.index = max(ghat$opt.index)
			ghat$opt.lambda = lambda.path[which(v==max(v[v<beta]))]
			ghat$network= LPGM.network(X,lambda=lambda.path,parallel=parallel,nCpus=nCpus,sym=sym,th=th)
			#ghat$network = lapply(1:length(ghat$network), function(i) WeightNet(ghat$network[[i]], ghat$D[[i]], thw=sth))
		}
		
		ghat$call <- match.call()
		cat(paste("\n", method, " Completed.", "\n", sep=""))
		class(ghat) <- "GMS"
		##return(ghat)
	}
	
	if(!is.null(ghat)){
		ghat$call <- match.call()
	}
	
	return(ghat)
}
