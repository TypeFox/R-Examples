# Kriging-based global sensitivity analysis taking into account both 
# the meta-model and the Monte-Carlo errors

# Author : loic le Gratiet, 2014

sobolGP <- function (
model,  
type="SK",
MCmethod="sobol",                                                                                                                                                               
X1,  
X2, 
nsim=100, 
nboot=1,
conf = 0.95,
sequential = FALSE, 
candidate = NULL, 
sequential.tot=FALSE,
max_iter = 1000
) 
{
    if(sequential){
    	ncandidate <- dim(candidate)[1]
    	dcandidate <- dim(candidate)[2]
    }	

    if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2))) 
        stop("The samples X1 and X2 must have the same dimensions")
    if(sequential){
	if(is.null(candidate)){
	  		stop("The number of candidate points must be greater than zero")
	}
    	if (ncol(X1) != ncol(candidate)){
	  		stop("The candidate points, X1 and X2 must have the same dimensions")
	}
    }
    if(MCmethod=="sobol"||MCmethod=="sobolEff"){
	if(sequential.tot){
	  stop("Sequential design for total indices is only available for sobol2002, sobol2007 and soboljansen methods")
	}
    }

    p <- ncol(X1)

    S <- list()

    if(sequential){
	Svar <- matrix(nrow = ncandidate, ncol = p)
    }
    if(sequential.tot){
	STotvar <- matrix(nrow = ncandidate, ncol = p)
    }

    output <- list()
    output$call$X1 <- X1
    output$call$X2 <- X2
    output$call$conf <- conf
    output$call$nboot <- nboot
    output$call$candidate <- candidate
    output$call$sequential <- sequential
    output$call$max_iter <- max_iter
    output$call$sequential.tot <- sequential.tot
    output$call$model <- model

Tot=FALSE
if(MCmethod=="sobol2002"||MCmethod=="sobol2007"||MCmethod=="soboljansen") Tot=TRUE

if(MCmethod!="sobol2007"){

    for (i in 1:p) {
        Xb <- X2
        Xb[, i] <- X1[, i]
        X <- rbind(X1, Xb)
	  nX <- dim(X)[1]

	if(sequential){
	  X <- rbind(X, data.frame(candidate))
	}
		
	rm(list=c("Xb"))

	  ysimu <- simulateGP.sobol(object = model, nsim = nsim,  newdata=X, 
                            cond=TRUE, checkNames=FALSE, max_iter=1000,type)

	if(MCmethod=="sobol"||MCmethod=="sobol2002"){
	  	S[[i]] <- sobolpickfreeze(ysimu[,1:(nX/2)] , ysimu[,(nX/2+1):nX],nboot)
	}
	if(MCmethod=="sobolEff"){
		S[[i]] <- sobolEffpickfreeze(ysimu[,1:(nX/2)] , ysimu[,(nX/2+1):nX],nboot)
	}
	if(MCmethod=="soboljansen"){
		S[[i]] <- soboljansenpickfreeze(ysimu[,1:(nX/2)] , ysimu[,(nX/2+1):nX],nboot)
	}

	  if(sequential){
	  	predCov <- predictGP.sobol(object = model, newdata1=data.frame(X), newdata2=data.frame(candidate), type=type, prednewdata1 = FALSE, prednewdata2 = TRUE)

	  	for(k in 1:ncandidate){
			ynew <- predCov$mean2[k]
			zsimu <- t(as.matrix(predCov$cov[-c((nX+1):(nX+ncandidate)),k])%*%rep(1,nsim))*(ynew-ysimu[,(nX+k)])/predCov$cov[(nX+k),k]+ysimu[,-c((nX+1):(nX+ncandidate))]

			if(MCmethod=="sobol"||MCmethod=="sobol2002"){
	  			Scand <- sobolpickfreeze(zsimu[,1:(nX/2)] , zsimu[,(nX/2+1):nX],nboot=1) 
			}
			if(MCmethod=="sobolEff"){
				Scand <- sobolEffpickfreeze(zsimu[,1:(nX/2)] , zsimu[,(nX/2+1):nX],nboot=1) 
			}
			if(MCmethod=="soboljansen"){
				Scand <- soboljansenpickfreeze(zsimu[,1:(nX/2)] , zsimu[,(nX/2+1):nX],nboot=1) 
			}

			rm(list=c("zsimu"))
			Svar[k,i] <- var(Scand)
			rm(list=c("Scand"))
	  	}
		rm(list=c("predCov", "ynew"))
	  }
	rm(list=c("ysimu" ))
    }

	

	if(Tot){
	
	    STot <- list()
	
    	for (i in 1:p) {
	        Xb <- X1
	        Xb[, i] <- X2[, i]
        	X <- rbind(X1, Xb)
		  nX <- dim(X)[1]
	  	if(sequential.tot){
		  X <- rbind(X, data.frame(candidate))
	  	}
			
		rm(list=c("Xb"))
	
		  ysimu <- simulateGP.sobol(object = model, nsim = nsim,  newdata=X, 
	                            cond=TRUE, checkNames=FALSE, max_iter=1000,type)
	
		if(MCmethod=="sobol2002"){ 
			STot[[i]] <- sobolT2002pickfreeze(ysimu[,1:(nX/2)],ysimu[,(nX/2+1):nX],nboot)
		}
		if(MCmethod=="soboljansen"){ 
			STot[[i]] <- sobolTjansenpickfreeze(ysimu[,1:(nX/2)],ysimu[,(nX/2+1):nX],nboot)
		}
	
		  if(sequential.tot){

		  	predCov <- predictGP.sobol(object = model, newdata1=data.frame(X), newdata2=data.frame(candidate), type=type)
	
		  	for(k in 1:ncandidate){
				ynew <- predCov$mean[k]
				zsimu <- t(as.matrix(predCov$cov[-c((nX+1):(nX+ncandidate)),k])%*%rep(1,nsim))*(ynew-ysimu[,(nX+k)])/predCov$cov[(nX+k),k]+ysimu[,-c((nX+1):(nX+ncandidate))]

				if(MCmethod=="sobol2002"){ 
					STotcand <- sobolT2002pickfreeze(zsimu[,1:(nX/2)],zsimu[,(nX/2+1):nX],nboot=1)
				}
				if(MCmethod=="soboljansen"){ 
					STotcand <- sobolTjansenpickfreeze(zsimu[,1:(nX/2)],zsimu[,(nX/2+1):nX],nboot=1)
				}

				rm(list=c("zsimu"))
				STotvar[k,i] <- var(STotcand)
				rm(list=c("STotcand"))
		  	}
		  rm(list=c("predCov", "ynew"))
	  	}
		rm(list=c("ysimu" ))
	    }
	}

		rm(list=c("X", "X1","X2"))
} 

if(MCmethod=="sobol2007"){
    STot <- list()

    for (i in 1:p) {
        Xb <- X1
        Xb[, i] <- X2[, i]
        X <- rbind(X1, Xb, X2)
	nX <- dim(X)[1]

	if(sequential||sequential.tot){
		X <- rbind(X, data.frame(candidate))
	}		

	rm(list=c("Xb"))

	  ysimu <- simulateGP.sobol(object = model, nsim = nsim,  newdata=X, 
                            cond=TRUE, checkNames=FALSE, max_iter=1000,type)
	

	S[[i]] <- sobol2007pickfreeze(ysimu[,1:(nX/3)] , ysimu[,(nX/3+1):(2*nX/3)] , ysimu[,(2*nX/3+1):nX],nboot)
	STot[[i]] <- sobolT2007pickfreeze(ysimu[,1:(nX/3)] ,ysimu[,(nX/3+1):(2*nX/3)] ,nboot)	

	  if(sequential||sequential.tot){

	  	predCov <- predictGP.sobol(object = model, newdata1=data.frame(X), newdata2=data.frame(candidate), type=type)

	  	for(k in 1:ncandidate){
			ynew <- predCov$mean[k]
			zsimu <- t(as.matrix(predCov$cov[-c((nX+1):(nX+ncandidate)),k])%*%rep(1,nsim))*(ynew-ysimu[,(nX+k)])/predCov$cov[(nX+k),k]+ysimu[,-c((nX+1):(nX+ncandidate))]

			if(sequential){
	  			Scand <- sobol2007pickfreeze(zsimu[,1:(nX/3)] , zsimu[,(nX/3+1):(2*nX/3)] , zsimu[,(2*nX/3+1):nX],nboot=1)
				Svar[k,i] <- var(Scand)
				rm(list=c("Scand"))
			}
			if(sequential.tot){
				STotcand <- sobolT2007pickfreeze(zsimu[,1:(nX/3)] ,zsimu[,(nX/3+1):(2*nX/3)], nboot=1 )	
				STotvar[k,i] <- var(STotcand)
				rm(list=c("STotcand"))
			}
			rm(list=c("zsimu"))

	  	}
		rm(list=c("predCov", "ynew"))
	  }
	rm(list=c("ysimu"))
    }
}

namesS <- c()
for (i in 1:p){
	namesS <- c(namesS,paste("S",i,sep=""))
}
names(S) <- namesS

if(Tot){
	namesStot <- c()
	for (i in 1:p){
		namesStot <- c(namesStot,paste("T",i,sep=""))
	}
	names(STot) <- namesStot
}


    output$S <- S
	rm(list=c("S"))

    if(sequential){
	   SumVar <- apply(Svar,1,sum)
	   output$S$xnew <- candidate[which.min(SumVar),]
	   output$S$xnewi <- which.min(SumVar)
		rm(list=c("Svar","SumVar"))
    }
   
    if(Tot){
	   output$T <- STot
		rm(list=c("STot"))
	   if(sequential.tot){
	   	SumVar <- apply(STotvar,1,sum)
	   	output$T$xnew <- candidate[which.min(SumVar),]
	   	output$T$xnewi <- which.min(SumVar)	
		rm(list=c("STotvar","SumVar"))	
	   }
    }

	class(output) <- "sobolGP"

	output$S$mean <- matrix(nrow = 1, ncol = p)
	output$S$var <- matrix(nrow = 1, ncol = p)
	output$S$ci <- matrix(nrow = 2, ncol = p)

	output$S$varPG <- matrix(nrow = 1, ncol = p)
	output$S$varMC <- matrix(nrow = 1, ncol = p)		

	if(Tot){
		output$T$mean <- matrix(nrow = 1, ncol = p)
		output$T$var <- matrix(nrow = 1, ncol = p)
		output$T$ci <- matrix(nrow = 2, ncol = p)

		output$T$varPG <- matrix(nrow = 1, ncol = p)
		output$T$varMC <- matrix(nrow = 1, ncol = p)	
	}
	for(i in 1:p){
		output$S$mean[1,i] <- mean(as.numeric(output$S[[i]]))
		output$S$var[1,i] <- var(as.numeric(output$S[[i]]))
		output$S$ci[1,i] <- quantile(as.numeric(output$S[[i]]), (1-conf)/2) 
		output$S$ci[2,i] <- quantile(as.numeric(output$S[[i]]), (1+conf)/2) 

		if(nboot==1){
			output$S$varPG[1,i] <- var(output$S[[i]])
		} else {
			output$S$varPG[1,i] <- mean(apply(output$S[[i]],1,var))
			output$S$varMC[1,i] <- mean(apply(output$S[[i]],2,var))
		}

		if(Tot){
			output$T$mean[1,i] <- mean(as.numeric(output$T[[i]]))
			output$T$var[1,i] <- var(as.numeric(output$T[[i]]))
			output$T$ci[1,i] <- quantile(as.numeric(output$T[[i]]), (1-conf)/2) 
			output$T$ci[2,i] <- quantile(as.numeric(output$T[[i]]), (1+conf)/2) 

			if(nboot==1){
				output$T$varPG[1,i] <- var(output$T[[i]])
			} else {
				output$T$varPG[1,i] <- mean(apply(output$T[[i]],1,var))
				output$T$varMC[1,i] <- mean(apply(output$T[[i]],2,var))
			}
		}
	}

	output$call$tot <- Tot
	output$call$method <- MCmethod
	output$call$type <- type
	output$call$nsim <- nsim
	

	return(output)
}



