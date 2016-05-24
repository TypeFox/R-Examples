glm.network <-
function(X, Y=NULL, link="gaussian", lambda, parallel=FALSE, nCpus = 4, standardize=TRUE,sym=TRUE,th=0){
	if(is.null(Y)){
		#-- cat("here 1\n")
		Z <- X
		p <- nrow(Z)
		q <- 0
	}
	
	if(!is.null(Y)){
		if(ncol(X) == ncol(Y)){
			Z <- rbind(X, Y)
			p = nrow(X)
			q = nrow(Y)
		}
	}
	

	if(length(lambda)>1){
		#ghat = array(0,dim=c(nrow(X),nrow(X),length(lambda)))
		ghat = array(0,dim=c(nrow(Z),nrow(Z),length(lambda)))
		#-- cat("here \n\n")
		wrapper1 <- function(i){
      #-- cat(i, "\n")
			tryCatch(
				{
					fit = glmnet::glmnet(t(Z[-i,]),Z[i,],family=link,lambda=lambda,standardize=standardize)
				},
				error = function(e) {
					fit = glmnetEmpty(t(Z[-i,]), lambda)
				})
			
			fit$beta=as.matrix(fit$beta)
			
			if(ncol(fit$beta)<length(lambda)){
				tmp = matrix(0,nrow = nrow(fit$beta),ncol = length(lambda))
				tmp[,1:ncol(fit$beta)]=fit$beta
				tmp[,ncol(fit$beta):length(lambda)] = fit$beta[,ncol(fit$beta)]
				fit$beta = tmp
			}

      
			if(i==1){
				ghat[i,2:nrow(Z),]=fit$beta
			}else if(i==nrow(Z)){
				ghat[i,1:(nrow(Z)-1),]=fit$beta
			}else{
				ghat[i,1:(i-1),]=fit$beta[1:(i-1),]
				ghat[i,(i+1):nrow(Z),]=fit$beta[i:nrow(fit$beta),]	
			}
	
			return(ghat[i,,])
		}

		if(parallel){
			#require('multicore')
			if(q == 0){
				ghat2 = parallel::mclapply(1:nrow(Z),wrapper1)
				#-- if(!sym){
					for(i in 1:nrow(Z)){   ghat[i,,]=ghat2[[i]] }
				#-- }
				if(sym){
					#-- for(i in 1:nrow(Z)){   ghat[i,,] = ANDNet(ghat2[[i]], th=th)  }
				  for(i in 1:length(lambda)){   ghat[,,i] = ANDNet(ghat[,,i], th=th)  }
				}
			}
			if(q != 0){
				ghat2 = parallel::mclapply((p+1):nrow(Z),wrapper1)
	# wooi: this might need to change ...because 1:p is nothing?!
				#-- if(!sym){
					for(i in (p+1):nrow(Z)){ ghat[i,,]=ghat2[[i-p]] }
				#-- }
				if(sym){
					#-- for(i in (p+1):nrow(Z)){ ghat[i,,]= ANDNet(ghat2[[i-p]], th=th)}
				  for(i in 1:length(lambda)){ ghat[,,i]= ANDNet(ghat[,,i], th=th)}
				}
			}
			return(ghat)
		}
	
		if(parallel==FALSE)
		{
			if(q == 0){
				ghat2=lapply(1:nrow(Z),wrapper1)	
				#-- if(!sym){
					for(i in 1:nrow(Z)){   ghat[i,,]=ghat2[[i]]  }
				#-- }
				if(sym){
					#-- for(i in 1:nrow(Z)){   ghat[i,,]= ANDNet(ghat2[[i]], th=th)  }
				  for(i in 1:length(lambda)){   ghat[,,i]= ANDNet(ghat[,,i], th=th)  }
				}
			}
			if(q != 0){
				ghat2=lapply((p+1):nrow(Z),wrapper1)
	# wooi: this might need to change ...because 1:p is nothing?!
				#-- if(!sym){
					for(i in  (p+1):nrow(Z)){  ghat[i,,]=ghat2[[i-p]]  }
				#-- }
				if(sym){
					#-- for(i in  (p+1):nrow(Z)){  ghat[i,,]= ANDNet(ghat2[[i-p]], th=th)  }
				  for(i in  1:length(lambda)){  ghat[,,i]= ANDNet(ghat[,,i], th=th)  }
				}
			}
			return(ghat)
		}
	}
	
	if(length(lambda) ==1){		
		ghat=matrix(0,nrow=nrow(Z),ncol=nrow(Z))
		
		if(parallel){
			###library(snowfall)
			snowfall::sfInit(parallel=TRUE, cpus=nCpus)
			
			snowfall::sfExport("X",local=TRUE)
			snowfall::sfExport("ghat",local=TRUE)
			snowfall::sfLibrary(glmnet::glmnet)
			
			wrapper2 <- function(i){
				tryCatch(
					{
						fit = glmnet::glmnet(t(Z[-i,]),Z[i,],family=link,lambda= lambda,standardize=standardize)
					},
					error = function(e) {
						fit = glmnetEmpty(t(Z[-i,]), lambda)
					}
				)

				fit$beta=as.numeric(fit$beta)
				if(i==1){
					ghat[i,2:nrow(Z)]=fit$beta
				}else if(i==nrow(Z)){
					ghat[i,1:(nrow(Z)-1)]=fit$beta
				}else{
					ghat[i,1:(i-1)]=fit$beta[1:(i-1)]
					ghat[i,(i+1):nrow(Z)]=c(fit$beta[i:length(fit$beta)])	
				}
				return(ghat[i,])
			}
			
			if(q == 0){
				snowfall::sfExport("wrapper2")
				ghat=snowfall::sfSapply(1:nrow(Z),wrapper2)
				snowfall::sfStop()				
			}
			if(q != 0){
				snowfall::sfExport("wrapper2")
				ghat=snowfall::sfSapply((p+1):nrow(Z),wrapper2)	
				snowfall::sfStop()	
			}
			if(sym){
				ghat = ANDNet(ghat, th=th)
			}
			return(ghat)
		}
		# wooi question: should run this again if parallel==T?
		if(parallel == FALSE){
			st = p+1
			if(q == 0){
				st = 1
			}
			for(i in st:nrow(Z)){
				#cat(i)
				tryCatch(
					{
						fit = glmnet::glmnet(t(Z[-i,]),Z[i,],family=link,lambda= lambda,standardize=standardize)
					},
					error = function(e) {
						fit = glmnetEmpty(t(Z[-i,]), lambda)
					}
				)
				fit$beta=as.numeric(fit$beta)
				if(i==1){
					ghat[i,2:nrow(Z)]=fit$beta
				}else if(i==nrow(Z)){
					ghat[i,1:(nrow(Z)-1)]=fit$beta
				}else{
					ghat[i,1:(i-1)]=fit$beta[1:(i-1)]
					ghat[i,(i+1):nrow(Z)]=c(fit$beta[i:length(fit$beta)])	
				}
			}
			if(sym){
				ghat = ANDNet(ghat, th=th)
			}
			return(ghat)
		}
	}
}
