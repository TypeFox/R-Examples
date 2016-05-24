glm.generic <-
function(X, method="GGM", link="gaussian", stability="bootstrap", N=100,beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL ,parallel=TRUE,nCpus=4,sym=TRUE,th=0,sth=0.8){
	##require('huge')
	##require('glmnet')
	
	if(is.null(lambda.path) ){
		#lmax = myglmnet.max(X)
		lmax = myglmnet.max(X, link=link)
 		lambda.path = exp(seq(log(lmax),log(lmin),l=nlams));
	}
	
	if(parallel == TRUE){
		if(stability == "bootstrap"){
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
			cat(paste(method, ": Conducting sampling ... in progress: ", floor(100*(i/N)), "%", collapse=""),"\r")
			flush.console()
			
			glmpois.good <- 1
			
			while(glmpois.good){
				# Make sure sample with no gene with all zero values
				good <- 1
				while(good){
					#-- index = sample(1:ncol(X),b,replace=F)
					index = sample(1:ncol(X),b,replace=replaceF)
					#-- if(sum(apply(X[,index], 1, sum)==0)==0){
					if(sum(apply(X[,index], 1, function(x) length(unique(x))==1))==0){
						good <- 0
					}
				}
				
				tryCatch(
						{
							#ghat.path$raw= glmpois(X[,index],lambda=lambda.path,parallel=T,nCpus=nCpus)
							ghat.path$raw= glm.network(X[,index], NULL, link=link, lambda=lambda.path,parallel=TRUE,nCpus=nCpus,sym=sym,th=th)	
							glmpois.good <- 0
						},
						error = function(e) {
							cat("glmnet returns empty model. Try again.")
						}
				)
			}
			
			for(j in 1:length(lambda.path)){
				tmp=ghat.path$raw[,,j]
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
			ghat$network = glm.network(X, NULL, link=link, lambda=lambda.path,parallel=TRUE,nCpus=nCpus,sym=sym,th=th)
			ghat$network = lapply(1:nlams,function(r){return(ghat$network[,,r])})
			ghat$network = lapply(1:ghat$network,function(r) WeightNet(ghat$network[[i]], ghat$D[[i]], thw=sth))		
		}
		else{
			ghat$opt.index = which(v==max(v[v<beta]))
			ghat$opt.index = max(ghat$opt.index)
			ghat$opt.lambda = lambda.path[which(v==max(v[v<beta]))]
			ghat$network = glm.network(X, NULL, link=link, lambda=lambda.path,parallel=TRUE,nCpus=nCpus,sym=sym,th=th)
			ghat$network =lapply(1:nlams,function(r){return(ghat$network[,,r])})		
		}
				
		ghat$call <- match.call()
		cat(paste("\n", method, " Completed.", "\n", sep=""))
		class(ghat) <- "GMS"
		return(ghat)
	}
	
	
	if(parallel == FALSE){
		
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
		ghat.path$path=vector("list",length(lambda.path))
		ghat.path$D = list()
		v=c()
		
		#-- ghat=list()
		#-- ghat.path$D = list()
		#-- v=c()
		
		for( j in 1:length(lambda.path)){
			#cat ("j=", j, " \t")
			cat(paste(method, ": Conducting sampling ... in progress: ", floor(100*(j/length(lambda.path))), "%", collapse=""),"\r")
			flush.console()
			D=matrix(0,nrow=nrow(X),ncol=nrow(X))
			
			for(i in 1:N){
				#cat("\n i=", i, "\t")
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
							#tmp=glmpois(X[,index],lambda=lambda.path[j],parallel=F)
							tmp=glm.network(X[,index],NULL, link=link, lambda=lambda.path[j],parallel=FALSE,sym=sym,th=th)
							glmpois.good <- 0
						},
						error = function(e) {
							cat("glmnet returns empty model. Try again.\n")
						}
					)
				} 
				
				tmp[abs(tmp)<1e-06]=0
				tmp[abs(tmp)>1e-06]=1
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
			
			#D=D/N
			#D=2*D*(1-D)
			#v=c(v,mean(D[upper.tri(D)]))			
		}
		
		v=cummax(v)
		ghat$v=v
		ghat$lambda.path = lambda.path
		ghat$D = ghat.path$D
		
		if(stability=="bootstrap"){
			ghat$opt.index = 1
			ghat$opt.lambda = lambda.path[1]
			ghat$network = glm.network(X, NULL, link=link, lambda=lambda.path,parallel=parallel,nCpus=nCpus,sym=sym,th=th)
			ghat$network = lapply(1:nlams,function(r){return(ghat$network[,,r])})
			ghat$network = lapply(1:ghat$network,function(r) WeightNet(ghat$network[[i]], ghat$D[[i]], thw=sth))	
		}
		else{
			ghat$opt.index = which(v==max(v[v<beta]))
			ghat$opt.index = max(ghat$opt.index)
			ghat$opt.lambda = lambda.path[which(v==max(v[v<beta]))]
			ghat$network = glm.network(X,NULL, link=link, lambda=lambda.path,parallel=parallel,nCpus=nCpus,sym=sym,th=th)
			ghat$network =lapply(1:nlams,function(r){return(ghat$network[,,r])})
		}
		
		ghat$call <- match.call()
		cat(paste("\n", method, " Completed.", "\n", sep=""))
		class(ghat) <- "GMS"
		return(ghat)
	}
		
}
