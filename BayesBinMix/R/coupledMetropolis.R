#library('label.switching')
#library('doParallel')
#library('foreach')

myDirichlet <- function(alpha){
	k <- length(alpha)
	theta <- rgamma(k, shape = alpha, rate = 1)
	return(theta/sum(theta))
}



toSolve <- function(a,n,p0){
	myVals <- exp(lgamma(2*a) - lgamma(a) + lgamma(a+n) - lgamma(2*a + n) + log(2))
	l <- abs(myVals - p0)
	return(a[which(l == min(l))][1])

}




complete.loglikelihood<-function(x,z,pars){
	g <- dim(pars)[1]
	n <- dim(x)[1]
	J <- dim(pars)[2]
	d <- dim(x)[2]
	logl <- rep(0, n)
	logpi <- log(pars[,J])
	z <- as.numeric(z)
	for (j in 1:d){
		logl <- logl + dbinom(x[,j], size = 1, prob = pars[,j][z], log = TRUE) 
	}
	logl <- logl + logpi[z]
	return(sum(logl))
}


gibbsBinMix <- function(alpha,beta,gamma,K,m,burn,data,thinning,z.true,outputDir){
	if (missing(K)) {stop(cat(paste("    [ERROR]: number of clusters (K) not provided."), "\n"))}
	if (missing(m)) {stop(cat(paste("    [ERROR]: number of MCMC iterations (m) not provided."), "\n"))}
	if (missing(thinning)) {thinning <- 5}
	if (K < 2) {
		stop(cat(paste("    [ERROR]: number of clusters (K) should be at least equal to 2."), "\n"))
	}
	d <- dim(data)[2]
	if (burn > m - 1) {
		stop(cat(paste("    [ERROR]: burn-in period (burn) not valid"), "\n"))
	}
	if (burn < 0) {
		stop(cat(paste("    [ERROR]: burn-in period (burn) not valid"), "\n"))
	}
	if (thinning > m) {
		stop(cat(paste("    [ERROR]: thinning not valid."), "\n"))
	}
	if (m < 1) {
		stop(cat(paste("    [ERROR]: mcmc iterations (m) not valid."), "\n"))
	}
	if (missing(alpha)) {alpha <- 1}
	if (missing(beta)) {beta <- 1}
	if (missing(gamma)) {gamma <- rep(1,K)}
	if (missing(data)) {stop(cat(paste("    [ERROR]: data is missing."), "\n"))}
	n <- dim(data)[1]
	x <- data

	dir.create(outputDir)
	setwd(outputDir)


	p <- myDirichlet(rep(1,K)) # initial values for cluster probabilities
	theta <- array(data = runif(K*d), dim = c(K,d)) # initial values for response probabilities
	z <- sample(K,n,replace=TRUE,prob = p)
	s <- l <- newL <- numeric(K)
	sx <- array(data = 0, dim = c(K,d))
	#output
	theta.file <- "theta.txt"
	z.file <- "z.txt"
	p.file <- "p.txt"
	conTheta = file(theta.file,open = "w")
	conZ = file(z.file,open = "w")
	conP = file(p.file,open = "w")

	for (iter in 1:m){
		#update cluster allocations (z)
		s <- rep(0,K)
		sx <- array(data = 0, dim = c(K,d))
		for(i in 1:n){
			for(k in 1:K){
				#l[k] <-  p[k]*prod((theta[k,]^x[i,])*((1-theta[k,])^(1-x[i,])))
				l[k] <-  log(p[k]) + sum(x[i,]*log(theta[k,]) + (1-x[i,])*log(1-theta[k,]))
			}
			for(k in 1:K){
				newL[k] <-  1/sum(exp(l - l[k]))
			}
			z[i] <- sample(K,1,prob = newL)
			#z[i] <- sample(K,1,prob = l)		
			s[z[i]] <- s[z[i]] + 1
			sx[z[i],] <- sx[z[i],] + x[i,]
		}
		p <- myDirichlet(gamma + s)
		for(k in 1:K){
			theta[k,] <- rbeta(d, shape1 = alpha + sx[k,], shape2 = beta + s[k] - sx[k,])
		}
		if((iter %% thinning == 0)&(iter > burn)){
			cat(theta,"\n",file=conTheta)
			cat(z,"\n",file=conZ)
			cat(p,"\n",file=conP)		
		}
		if(iter %% (m/100) == 0){
			cat(paste("Running Gibbs: ",100*round(iter/m,3),"% completed.",sep=""),"\n");
		}
	}
	close(conTheta)
	close(conP)
	close(conZ)
	cat(paste("Gibbs finished."),"\n")
	cat("\n")
	cat(paste("Dealing with Label Switching..."),"\n")
	tt <- read.table("theta.txt")
	m <- dim(tt)[1]
	d <- dim(x)[2]
	J <- d + 1
	mcmc <- array(data = NA, dim = c(m,K,J))
	for (j in 1:d){
		for(k in 1:K){
			mcmc[,k,j] <- tt[,(j-1)*K + k]
		}
	}
	tt <- read.table("p.txt")
	for(k in 1:K){
		mcmc[,k,J] <- tt[,k]
	}
	allocations <- as.matrix(read.table("z.txt",as.is = TRUE))
	m <- dim(mcmc)[1]
	iter <- 1
	ll <- complete.loglikelihood(x,allocations[iter,],mcmc[iter,,])
	maxLL <- ll
	maxIter <- iter
	for (iter in 1:m){
		ll <- complete.loglikelihood(x,allocations[iter,],mcmc[iter,,])
		if(ll > maxLL){maxLL <- ll;maxIter <- iter;cat(paste("Found new Complete MLE: ", ll,sep=""),"\n")}
	}
	# classification probs
	pMatrix <- array(data = NA, dim = c(m,n,K))
	for (iter in 1:m){
		p <- mcmc[iter,,J]
		theta <- array( mcmc[iter,,1:(J-1)],dim = c(K,J-1) )
		for(i in 1:n){
			for(k in 1:K){
				#pMatrix[iter,i,k] <-  p[k]*prod((theta[k,]^x[i,])*((1-theta[k,])^(1-x[i,])))
				l[k] <- log(p[k]) + sum(x[i,]*log(theta[k,]) + (1-x[i,])*log(1-theta[k,]))
			}
			#pMatrix[iter,i,] <- pMatrix[iter,i,]/sum(pMatrix[iter,i,])
			for(k in 1:K){
				pMatrix[iter,i,k] <- 1/sum(exp(l - l[k]))
				if(is.na(pMatrix[iter,i,k]) == TRUE){pMatrix[iter,i,k] = 0}
			}
		}
		if(iter %% 1000 == 0){cat(paste(" classification probs: ",100*round(iter/m,3),"% completed",sep=""),"\n");}
	}

	if(missing(z.true)==TRUE){
		ls <- label.switching( K = K, method = c("STEPHENS","ECR","ECR-ITERATIVE-1"),
					zpivot = allocations[maxIter,], z = allocations, complete = complete.loglikelihood, data = x,
					mcmc = mcmc,prapivot = mcmc[maxIter,,], p = pMatrix)
	}else{
		ls <- label.switching( K = K, method = c("STEPHENS","ECR","ECR-ITERATIVE-1"),
					zpivot = allocations[maxIter,], z = allocations, complete = complete.loglikelihood, data = x,
					mcmc = mcmc,p = pMatrix,groundTruth = z.true)
	}
	#reordering allocations
	allocationsECR <- allocationsKL <- allocationsECR.ITERATIVE1 <- allocations
	for (i in 1:m){
		myPerm <- order(ls$permutations$"ECR"[i,])
		allocationsECR[i,] <- myPerm[allocations[i,]]
		myPerm <- order(ls$permutations$"STEPHENS"[i,])
		allocationsKL[i,] <- myPerm[allocations[i,]]
		myPerm <- order(ls$permutations$"ECR-ITERATIVE-1"[i,])
		allocationsECR.ITERATIVE1[i,] <- myPerm[allocations[i,]]

	}
	write.table(allocationsECR, file = "z.ECR.txt")
	write.table(allocationsKL, file = "z.KL.txt")
	write.table(allocationsECR.ITERATIVE1, file = "z.ECR-ITERATIVE1.txt")



	reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"ECR")$output
	tt <- read.table("theta.txt")
	write.table(reordered.mcmc, file = "reorderedMCMC-ECR.txt",col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
	reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"ECR-ITERATIVE-1")$output
	write.table(reordered.mcmc, file = "reorderedMCMC-ECR-ITERATIVE1.txt",col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
	reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"STEPHENS")$output
	write.table(reordered.mcmc, file = "reorderedMCMC-STEPHENS.txt",col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
	write.table(mcmc, file = "rawMCMC.txt",col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
	write.table(file = "reorderedSingleBestClusterings.txt",t(ls$clusters[c(1,2,3),]),row.names = paste("z",1:n,sep="."))
	file.remove("p.txt")
	file.remove("theta.txt")
	cat(paste("raw MCMC parameters written to: \'rawMCMC.txt\' "),"\n")
	cat(paste("raw MCMC latent allocations written to: \'z.txt\' "),"\n")
	cat(paste("reordered MCMC output written to: "),"\n")
	cat(paste("     (Method 1):     \'reorderedMCMC-ECR.txt\'"),"\n")
	cat(paste("     (Method 2):     \'reorderedMCMC-ECR-ITERATIVE1.txt\'"),"\n")
	cat(paste("     (Method 3):     \'reorderedMCMC-STEPHENS.txt\'"),"\n")
	cat(paste("reordered single best clusterings written to: \'reorderedSingleBestClusterings.txt\' "),"\n")
	cat(paste("reordered MCMC latent allocations written to: "),"\n")
	cat(paste("     (Method 1):     \'z.ECR.txt\'"),"\n")
	cat(paste("     (Method 2):     \'z.KL.txt\'"),"\n")
	cat(paste("     (Method 3):     \'z.ECR-ITERATIVE1.txt\'"),"\n")
	setwd("../")

}



collapsedGibbsBinMix <- function(alpha,beta,gamma,K,m,burn,data,thinning,z.true,outputDir){
	if (missing(K)) {stop(cat(paste("    [ERROR]: number of clusters (K) not provided."), "\n"))}
	if (missing(m)) {stop(cat(paste("    [ERROR]: number of MCMC iterations (m) not provided."), "\n"))}
	if (missing(thinning)) {thinning <- 5}
	d <- dim(data)[2]
	if (burn > m - 1) {
		stop(cat(paste("    [ERROR]: burn-in period (burn) not valid"), "\n"))
	}
	if (burn < 0) {
		stop(cat(paste("    [ERROR]: burn-in period (burn) not valid"), "\n"))
	}
	if (thinning > m) {
		stop(cat(paste("    [ERROR]: thinning not valid."), "\n"))
	}
	if (m < 1) {
		stop(cat(paste("    [ERROR]: mcmc iterations (m) not valid."), "\n"))
	}
	if (missing(alpha)) {alpha <- 1}
	if (missing(beta)) {beta <- 1}
	if (missing(gamma)) {gamma <- rep(1,K)}
	if (missing(data)) {stop(cat(paste("    [ERROR]: data is missing."), "\n"))}
	d <- dim(data)[2]
	n <- dim(data)[1]
	x <- data

	dir.create(outputDir)
	setwd(outputDir)

	p <- myDirichlet(rep(1,K)) # initial values for cluster probabilities
	theta <- array(data = runif(K*d), dim = c(K,d)) # initial values for response probabilities
	z <- sample(K,n,replace=TRUE,prob = p)

	s <- l <- newL <- numeric(K)
	sx <- array(data = 0, dim = c(K,d))
	#output
	theta.file <- "theta.collapsed.txt"
	z.file <- "z.collapsed.txt"
	p.file <- "p.collapsed.txt"
	conTheta = file(theta.file,open = "w")
	conZ = file(z.file,open = "w")
	conP = file(p.file,open = "w")
#	initialization with standard Gibbs 
	for (iter in 1:100){
		s <- rep(0,K)
		sx <- array(data = 0, dim = c(K,d))
		for(i in 1:n){
			for(k in 1:K){
				l[k] <-  log(p[k]) + sum(x[i,]*log(theta[k,]) + (1-x[i,])*log(1-theta[k,]))
			}
			for(k in 1:K){
				newL[k] <-  1/sum(exp(l - l[k]))
			}
			z[i] <- sample(K,1,prob = newL)
			s[z[i]] <- s[z[i]] + 1
			sx[z[i],] <- sx[z[i],] + x[i,]
		}
		p <- myDirichlet(gamma + s)
		for(k in 1:K){
			theta[k,] <- rbeta(d, shape1 = alpha + sx[k,], shape2 = beta + s[k] - sx[k,])
		}
	}
#	collapsed sampler
	nMatrix <- array(data = 0, dim = c(n,K))
	sMatrix <- array(data = 0, dim = c(n,K,d))

	for (i in 1:n){
		nMatrix[i,] <- s
		sMatrix[i,,] <- sx
	}
	for (i in 1:n){
		nMatrix[i,z[i]] <- nMatrix[i,z[i]] - 1
		for (j in 1:d){
			sMatrix[i,z[i],j] <- sMatrix[i,z[i],j] - x[i,j]
		}
	}

	reallocationAcceptanceRatio <- 0
	reallocationAcceptanceRatio2 <- 0
	a1 <- a2 <- 1 #parameters of Beta for the reallocation proposal
	lar <- c()
	for(iter in 1:m){
		zOld <- z
		for(i in 1:n){
			nMatrix[i,1:K] <- s[1:K]
			sMatrix[i,1:K,] <- sx[1:K,]
			nMatrix[i,z[i]] <- s[z[i]] - 1
			sMatrix[i,z[i],] <- sx[z[i],] - x[i,]

			A1 <- which(x[i,] == 1)
			A0 <- which(x[i,] == 0)
			for(k in 1:K){
				newL[k] <- log(nMatrix[i,k] + gamma[k]) - d*log(alpha + beta + nMatrix[i,k])
				if(length(A1 > 0)){
					newL[k] <- newL[k] + sum(log(alpha + sMatrix[i,k,A1]))
				}
				if(length(A0 > 0)){
					newL[k] <- newL[k] + sum(log(beta + nMatrix[i,k] - sMatrix[i,k,A0]))
				}
			}
			newL[1:K] <- exp(newL[1:K])
			z[i] <- sample(K,1,prob = newL[1:K])
			#if( zOld[i] != z[i] ){
			sx[zOld[i],] <- sx[zOld[i],] - x[i,]
			sx[z[i],] <- sx[z[i],] + x[i,]
			s[zOld[i]] <- s[zOld[i]] - 1
			s[z[i]] <- s[z[i]] + 1
		}


#		reallocation proposal 1
		if (K > 1){
			myPair <- sample(K,2,replace = FALSE)
			propZ <- z
			set1 <- which(z == myPair[1])
			set2 <- which(z == myPair[2])
			myP <- rbeta(1,shape1 = gamma[myPair[1]],shape2 = gamma[myPair[2]])
			propZ[c(set1,set2)] <- myPair[sample(2,length(set1) + length(set2), replace = TRUE,prob = c(myP,1-myP) )]
			newSet1 <- which(propZ == myPair[1])
			newSet2 <- which(propZ == myPair[2])
			nNew <- c(length(newSet1),length(newSet2))
			nOld <- c(length(set1),length(set2))
			if(nOld[1] == 0){sOld1 <- rep(0,d)}else{sOld1 <- colSums(array(x[set1,],dim = c(nOld[1],d)))}
			if(nOld[2] == 0){sOld2 <- rep(0,d)}else{sOld2 <- colSums(array(x[set2,],dim = c(nOld[2],d)))}
			sOld <- rbind(sOld1,sOld2)
			if(nNew[1] == 0){sNew1 <- rep(0,d)}else{sNew1 <- colSums(array(x[newSet1,],dim = c(nNew[1],d)))}
			if(nNew[2] == 0){sNew2 <- rep(0,d)}else{sNew2 <- colSums(array(x[newSet2,],dim = c(nNew[2],d)))}
			sNew <- rbind(sNew1,sNew2)
			logAR <- 0
			for (i in 1:2){
				logAR <- logAR + d*(lgamma(alpha + beta + nOld[i]) - lgamma(alpha + beta + nNew[i]) ) + sum(lgamma(alpha + sNew[i,]) + lgamma(beta + nNew[i] - sNew[i,]) - lgamma(alpha + sOld[i,]) - lgamma(beta + nOld[i] - sOld[i,]))
			}
			if( log(runif(1)) < logAR ){
				reallocationAcceptanceRatio <- reallocationAcceptanceRatio + 1 
				z <- propZ
				s <- rep(0,K)
				sx <- array(data = 0, dim = c(K,d))
				for(i in 1:n){
					s[z[i]] <- s[z[i]] + 1
					sx[z[i],] <- sx[z[i],] + x[i,]
				}
				for(i in 1:n){
					nMatrix[i,] <- s
					sMatrix[i,,] <- sx
				}
				for (i in 1:n){
					nMatrix[i,z[i]] <- nMatrix[i,z[i]] - 1
					for (j in 1:d){
						sMatrix[i,z[i],j] <- sMatrix[i,z[i],j] - x[i,j]
					}
				}
			}
		}
#		end of reallocation proposal 1

#		reallocation proposal2
		if(K > 1){
			myPair <- sample(K,2,replace = FALSE)
			set1 <- which(z == myPair[1])
			set2 <- which(z == myPair[2])
			nOld <- c(length(set1),length(set2))
			propZ <- z
			if(nOld[1] > 0){
				randomSize <- 1 + floor(nOld[1]*runif(1))
				randomIndex <- set1[sample(nOld[1],randomSize,replace = FALSE)]
				propZ[randomIndex] <- rep(myPair[2],randomSize)
				newSet1 <- which(propZ == myPair[1])
				newSet2 <- which(propZ == myPair[2])
				nNew <- c(length(newSet1),length(newSet2))
				if(nOld[1] == 0){sOld1 <- rep(0,d)}else{sOld1 <- colSums(array(x[set1,],dim = c(nOld[1],d)))}
				if(nOld[2] == 0){sOld2 <- rep(0,d)}else{sOld2 <- colSums(array(x[set2,],dim = c(nOld[2],d)))}
				sOld <- rbind(sOld1,sOld2)
				if(nNew[1] == 0){sNew1 <- rep(0,d)}else{sNew1 <- colSums(array(x[newSet1,],dim = c(nNew[1],d)))}
				if(nNew[2] == 0){sNew2 <- rep(0,d)}else{sNew2 <- colSums(array(x[newSet2,],dim = c(nNew[2],d)))}
				sNew <- rbind(sNew1,sNew2)
				logAR <- log(nOld[1]) - log(nOld[2] + randomSize) -( lgamma(nOld[1] - randomSize + 1) + lgamma(nOld[2]+ randomSize + 1) - lgamma(nOld[2] + 1) - lgamma(nOld[1] + 1))
				for(i in 1:2){
					logAR <- logAR + d*(lgamma(alpha + beta + nOld[i]) - lgamma(alpha + beta + nNew[i]) ) + sum(lgamma(alpha + sNew[i,]) + lgamma(beta + nNew[i] - sNew[i,]) - lgamma(alpha + sOld[i,]) - lgamma(beta + nOld[i] - sOld[i,]))
				}
				logAR <- logAR + sum(lgamma(gamma[myPair] + nNew)) - sum(lgamma(gamma[myPair] + nOld)) 
				if( log(runif(1)) < logAR ){
					reallocationAcceptanceRatio2 <- reallocationAcceptanceRatio2 + 1 
					z <- propZ
					s <- rep(0,K)
					sx <- array(data = 0, dim = c(K,d))
					for(i in 1:n){
						s[z[i]] <- s[z[i]] + 1
						sx[z[i],] <- sx[z[i],] + x[i,]
					}
					for(i in 1:n){
						nMatrix[i,] <- s
						sMatrix[i,,] <- sx
					}
					for (i in 1:n){
						nMatrix[i,z[i]] <- nMatrix[i,z[i]] - 1
						for (j in 1:d){
							sMatrix[i,z[i],j] <- sMatrix[i,z[i],j] - x[i,j]
						}
					}
				}
			}
		}

#		end of reallocation proposal 2
		if((iter %% thinning == 0)&(iter > burn)){
			for(k in 1:K){
				myIndex <- which(z == k)
				s[k] <- length(myIndex)
				if(s[k] > 1){
					sx[k,] <- colSums(x[myIndex,])
				}else{
					if(s[k] == 1){
						sx[k,] <- x[myIndex,]
					}else{
						sx[k,] <- rep(0,d)
					}

				}
				theta[k,] <- rbeta(d, shape1 = alpha + sx[k,], shape2 = beta + s[k] - sx[k,])
			}
			p <- myDirichlet(gamma + s)
		
			cat(theta,"\n",file=conTheta)
			cat(z,"\n",file=conZ)
			cat(p,"\n",file=conP)		
		}
		if(iter %% (m/100) == 0){cat(paste("Running collapsed Gibbs: ",100*round(iter/m,3),"% completed. Reallocation proposal acceptance rates: ", 100*round(reallocationAcceptanceRatio/iter,3),"%, ",100*round(reallocationAcceptanceRatio2/iter,3),"%.",sep=""),"\n");
		#matplot(read.table("p.collapsed.txt"),type = "l",lty = 1)
	}


	}

	close(conTheta)
	close(conP)
	close(conZ)
	cat(paste("Collapsed Gibbs finished."),"\n")
	cat("\n")
	if(K > 1){
		cat(paste("Dealing with Label Switching..."),"\n")
		tt <- read.table("theta.collapsed.txt")
		m <- dim(tt)[1]
		d <- dim(x)[2]
		J <- d + 1
		mcmc <- array(data = NA, dim = c(m,K,J))
		for (j in 1:d){
			for(k in 1:K){
				mcmc[,k,j] <- tt[,(j-1)*K + k]
			}
		}
		tt <- read.table("p.collapsed.txt")
		for(k in 1:K){
			mcmc[,k,J] <- tt[,k]
		}
		allocations <- as.matrix(read.table("z.collapsed.txt",as.is = TRUE))
		m <- dim(mcmc)[1]
		iter <- 1
		ll <- complete.loglikelihood(x,allocations[iter,],mcmc[iter,,])
		maxLL <- ll
		maxIter <- iter
		for (iter in 1:m){
			ll <- complete.loglikelihood(x,allocations[iter,],mcmc[iter,,])
			if(ll > maxLL){maxLL <- ll;maxIter <- iter;cat(paste("Found new Complete MLE: ", ll,sep=""),"\n")}
		}
		# classification probs
		pMatrix <- array(data = NA, dim = c(m,n,K))
		for (iter in 1:m){
			p <- mcmc[iter,,J]
			theta <- array( mcmc[iter,,1:(J-1)],dim = c(K,J-1) )
			for(i in 1:n){
				for(k in 1:K){
					l[k] <- log(p[k]) + sum(x[i,]*log(theta[k,]) + (1-x[i,])*log(1-theta[k,]))
				}
				for(k in 1:K){
					pMatrix[iter,i,k] <- 1/sum(exp(l - l[k]))
					if(is.na(pMatrix[iter,i,k]) == TRUE){pMatrix[iter,i,k] = 0}
				}
			}
			if(iter %% 1000 == 0){cat(paste(" classification probs: ",100*round(iter/m,3),"% completed",sep=""),"\n");}
		}


		if(missing(z.true)==TRUE){
			ls <- label.switching( method = c("STEPHENS","ECR","ECR-ITERATIVE-1"),
						zpivot = allocations[maxIter,], z = allocations,K = K, complete = complete.loglikelihood, data = x,
						prapivot = mcmc[maxIter,,], mcmc = mcmc, p = pMatrix)
		}else{
			ls <- label.switching( method = c("STEPHENS","ECR","ECR-ITERATIVE-1"),
						zpivot = allocations[maxIter,], z = allocations,K = K, complete = complete.loglikelihood, data = x,
						prapivot = mcmc[maxIter,,], mcmc = mcmc, p = pMatrix,groundTruth = z.true)
		}
		reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"ECR")$output

		tt <- read.table("theta.collapsed.txt")

		write.table(reordered.mcmc, file = "reorderedMCMC-ECR.collapsed.txt",col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
		reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"ECR-ITERATIVE-1")$output
		write.table(reordered.mcmc, file = "reorderedMCMC-ECR-ITERATIVE1.collapsed.txt",col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
		reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"STEPHENS")$output
		write.table(reordered.mcmc, file = "reorderedMCMC-STEPHENS.collapsed.txt",col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
		write.table(mcmc, file = "rawMCMC.collapsed.txt",col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
		write.table(file = "reorderedSingleBestClusterings.collapsed.txt",t(ls$clusters[c(1,2,3),]),row.names = paste("z",1:n,sep="."))
		file.remove("p.collapsed.txt")
		file.remove("theta.collapsed.txt")
		cat(paste("raw MCMC parameters written to: \'rawMCMC.collapsed.txt\' "),"\n")
		cat(paste("raw MCMC latent allocations written to: \'z.collapsed.txt\' "),"\n")
		cat(paste("reordered MCMC output written to: "),"\n")
		cat(paste("     (Method 1):     \'reorderedMCMC-ECR.collapsed.txt\'"),"\n")
		cat(paste("     (Method 2):     \'reorderedMCMC-ECR-ITERATIVE1.collapsed.txt\'"),"\n")
		cat(paste("     (Method 3):     \'reorderedMCMC-STEPHENS.collapsed.txt\'"),"\n")
		cat(paste("reordered single best clusterings written to: \'reorderedSingleBestClusterings.collapsed.txt\' "),"\n")
	}
	setwd("../")
}


######################################################




allocationSamplerBinMix <- function(Kmax, alpha,beta,gamma,m,burn,data,thinning,z.true,ClusterPrior,ejectionAlpha,Kstart,outputDir,metropolisMoves,reorderModels,heat,zStart,LS){
#	if(dir.exists(outputDir) == TRUE){
#		stop(cat(paste("    [ERROR]: directory exists, please provide different name."), "\n"))
#	}
	 #dir.create(outputDir)
	setwd(outputDir)
	cat(paste0("changing working directory to ",outputDir),"\n")
	if (missing(m)) {stop(cat(paste("    [ERROR]: number of MCMC iterations (m) not provided."), "\n"))}
	if (missing(thinning)) {thinning <- 5}
	if (missing(LS)) {LS <- TRUE}
	if (missing(thinning)) {ejectionAlpha <- 1}
	if (burn > m - 1) {
		stop(cat(paste("    [ERROR]: burn-in period (burn) not valid"), "\n"))
	}
	if (burn < 0) {
		stop(cat(paste("    [ERROR]: burn-in period (burn) not valid"), "\n"))
	}
	if (thinning > m) {
		stop(cat(paste("    [ERROR]: thinning not valid."), "\n"))
	}
	if (m < 1) {
		stop(cat(paste("    [ERROR]: mcmc iterations (m) not valid."), "\n"))
	}
	if (missing(alpha)) {alpha <- 1}
	if (missing(beta)) {beta <- 1}
	if (missing(data)) {stop(cat(paste("    [ERROR]: data is missing."), "\n"))}
	x <- data
	d <- dim(data)[2]
	n <- dim(data)[1]
	if (missing(Kmax)) {Kmax <- floor((d + 1)/2)}
	if (missing(gamma)) {gamma <- rep(1,Kmax)}
	if (length(table(gamma)) > 1){
		stop(cat(paste("    [ERROR]: Dirichlet prior parameters should be the same."), "\n"))
	}
	birthProbs <- rep(1,Kmax)
	birthProbs[2:(Kmax - 1)] <- rep(0.5,Kmax - 2)
	birthProbs[Kmax] <- 0
	deathProbs <- 1 - birthProbs

	# define uniform or truncated poisson(1) prior on K
	priorK <- numeric(Kmax)
	if (missing(ClusterPrior)) {stop(cat(paste("    [ERROR]: ClusterPrior not defined (uniform of poisson)."), "\n"))}
	if (missing(ejectionAlpha)) {stop(cat(paste("    [ERROR]: ejectionAlpha not defined (0 < a < 1)."), "\n"))}
	if (ejectionAlpha < 0) {stop(cat(paste("    [ERROR]: ejectionAlpha error (0 < a < 1)."), "\n"))}
	if (ejectionAlpha > 1) {stop(cat(paste("    [ERROR]: ejectionAlpha error (0 < a < 1)."), "\n"))}
	if(ClusterPrior == "uniform"){
		priorK <- rep(log(1/Kmax),Kmax)
	}
	if(ClusterPrior == "poisson"){
		denom <- log(ppois(Kmax, lambda = 1, lower.tail = TRUE) - dpois(0,lambda = 1))
		for (k in 1:Kmax){
			priorK[k] <- dpois(k,lambda = 1, log = TRUE) - denom 
		}
	}
	if(missing(Kstart)){
		Kstart <- 1
	}
	K <- Kstart
	p <- myDirichlet(rep(1,K)) # initial values for cluster probabilities
	theta <- array(data = runif(Kmax*d), dim = c(Kmax,d)) # initial values for response probabilities
	s <- l <- newL <- numeric(Kmax)
	sx <- array(data = 0, dim = c(Kmax,d))


	if(missing(zStart) == TRUE){
		cat(paste("Initializing..."),"\n")
		cat("\n")
		z <- sample(K,n,replace=TRUE,prob = p)
		cat(paste("Allocation sampler running..."),"\n")
		cat("\n")
	}else{
		z <- zStart
		s[1:K] <- rep(0,K)
		sx[1:K,] <- array(data = 0, dim = c(K,d))
		for(i in 1:n){
			s[z[i]] <- s[z[i]] + 1
			sx[z[i],] <- sx[z[i],] + x[i,]
		}

	}

	#output
	k.file <- "K.txt"
	theta.file <- "theta.varK.txt"
	z.file <- "z.varK.txt"
	p.file <- "p.varK.txt"
	conK = file(k.file,open = "w")
	conTheta = file(theta.file,open = "w")
	conZ = file(z.file,open = "w")
	conP = file(p.file,open = "w")
#	initialization with standard Gibbs 
	#print("here1")
	#par(mfrow = c(1,1))
	if(missing(zStart) == TRUE){
	for (iter in 1:5){
		s[1:K] <- rep(0,K)
		sx[1:K,] <- array(data = 0, dim = c(K,d))
		for(i in 1:n){
			for(k in 1:K){
				l[k] <-  log(p[k]) + sum(x[i,]*log(theta[k,]) + (1-x[i,])*log(1-theta[k,]))
			}
			for(k in 1:K){
				newL[k] <-  1/sum(exp(l[1:K] - l[k]))
			}
			if(is.nan(sum(newL[1:K]))==TRUE){newL[1:K]<-rep(1,K)}
			z[i] <- sample(K,1,prob = newL[1:K])
			s[z[i]] <- s[z[i]] + 1
			sx[z[i],] <- sx[z[i],] + x[i,]
		}
		p <- myDirichlet(heat*gamma[1:K] + heat*s[1:K] + 1 - heat)
		for(k in 1:K){
			theta[k,] <- rbeta(d, shape1 = heat*alpha + heat*sx[k,] + 1 - heat, shape2 = heat*beta + heat*s[k] - heat*sx[k,] + 1 - heat)
		}

		
	}
	}

#	collapsed sampler
	nMatrix <- array(data = 0, dim = c(n,Kmax))
	sMatrix <- array(data = 0, dim = c(n,Kmax,d))
	#print("here2")

	for (i in 1:n){
		nMatrix[i,1:K] <- s[1:K]
		sMatrix[i,1:K,] <- sx[1:K,]
	}
	for (i in 1:n){
		nMatrix[i,z[i]] <- nMatrix[i,z[i]] - 1
		for (j in 1:d){
			sMatrix[i,z[i],j] <- sMatrix[i,z[i],j] - x[i,j]
		}
	}


	
	reallocationAcceptanceRatio <- 0
	reallocationAcceptanceRatio2 <- 0
	reallocationAcceptanceRatio3 <- 0
	reallocationAcceptanceRatio4 <- 0
	a1 <- a2 <- 1 #parameters of Beta for the reallocation proposal
	pZERO <- ejectionAlpha
	constLGAMMA <- 2*lgamma(ejectionAlpha) - lgamma(2*ejectionAlpha)
	kValues <- c()
	#plot(c(1,m),c(0,1),type = "n")
	if(missing(metropolisMoves) == TRUE){metropolisMoves <- c('M1','M2','M3','M4')}
	lmm <- length(metropolisMoves)

	cat(paste("    Reallocation proposal acceptance rates: Move 1, Move 2, Move 3, Move 4"),"\n")
		                                              

	for(iter in 1:m){
		zOld <- z
		sOld <- s
		sxOld <- sx
		if(K > 1){
			myPair <- sample(K,2,replace = FALSE)
			j1 <- myPair[1]
			j2 <- myPair[2]
			set1 <- which(z == j1)
			set2 <- which(z == j2)
			zOld[set1] = rep(j2,length(set1))
			zOld[set2] = rep(j1,length(set2))
			z <- zOld
			sOld[j1] <- s[j2]
			sOld[j2] <- s[j1]
			sxOld[j1,] <- sx[j2,]
			sxOld[j2,] <- sx[j1,]
			s <- sOld
			sx <- sxOld
		}
		for(i in 1:n){
			nMatrix[i,1:K] <- s[1:K]
			sMatrix[i,1:K,] <- sx[1:K,]
			nMatrix[i,z[i]] <- s[z[i]] - 1
			sMatrix[i,z[i],] <- sx[z[i],] - x[i,]

			A1 <- which(x[i,] == 1)
			A0 <- which(x[i,] == 0)
			for(k in 1:K){
				newL[k] <- heat*log(nMatrix[i,k] + gamma[k]) - heat*d*log(alpha + beta + nMatrix[i,k])
				if(length(A1 > 0)){
					newL[k] <- newL[k] + heat*sum(log(alpha + sMatrix[i,k,A1]))
				}
				if(length(A0 > 0)){
					newL[k] <- newL[k] + heat*sum(log(beta + nMatrix[i,k] - sMatrix[i,k,A0]))
				}
			}
			newL[1:K] <- exp(newL[1:K])
			if( max(newL[1:K]) == 0){newL[1:K] <- rep(1,K);cat("oops","\n")}
			z[i] <- sample(K,1,prob = newL[1:K])
			#if( zOld[i] != z[i] ){
			sx[zOld[i],] <- sx[zOld[i],] - x[i,]
			sx[z[i],] <- sx[z[i],] + x[i,]
			s[zOld[i]] <- s[zOld[i]] - 1
			s[z[i]] <- s[z[i]] + 1
			#nMatrix[,z[i]] <- s[z[i]] - 1
			#sMatrix[,z[i],] <- sx[z[i],] - x[i,]
				#nMatrix[,zOld[i]] <- nMatrix[,zOld[i]] - 1
				#nMatrix[,z[i]] <- nMatrix[,z[i]] + 1
				#myMat <- matrix(rep(x[i,],n),nrow = n,byrow=TRUE)
				#sMatrix[,zOld[i],] <- sMatrix[,zOld[i],] - myMat
				#sMatrix[,z[i],] <- sMatrix[,z[i],] + myMat
			#}
		}
	mNumber <- sample(metropolisMoves,1)
	#mNumber <- 5


		if(mNumber == 'M1'){
#		if(1 > 2){
#		reallocation proposal 1
		if (K > 1){
			myPair <- sample(K,2,replace = FALSE)
			propZ <- z
			set1 <- which(z == myPair[1])
			set2 <- which(z == myPair[2])
			myP <- rbeta(1,shape1 = gamma[myPair[1]],shape2 = gamma[myPair[2]])
			propZ[c(set1,set2)] <- myPair[sample(2,length(set1) + length(set2), replace = TRUE,prob = c(myP,1-myP) )]
			newSet1 <- which(propZ == myPair[1])
			newSet2 <- which(propZ == myPair[2])
			nNew <- c(length(newSet1),length(newSet2))
			nOld <- c(length(set1),length(set2))
			if(nOld[1] == 0){sOld1 <- rep(0,d)}else{sOld1 <- colSums(array(x[set1,],dim = c(nOld[1],d)))}
			if(nOld[2] == 0){sOld2 <- rep(0,d)}else{sOld2 <- colSums(array(x[set2,],dim = c(nOld[2],d)))}
			sOld <- rbind(sOld1,sOld2)
			if(nNew[1] == 0){sNew1 <- rep(0,d)}else{sNew1 <- colSums(array(x[newSet1,],dim = c(nNew[1],d)))}
			if(nNew[2] == 0){sNew2 <- rep(0,d)}else{sNew2 <- colSums(array(x[newSet2,],dim = c(nNew[2],d)))}
			sNew <- rbind(sNew1,sNew2)
			logAR <- 0
			for (i in 1:2){
				logAR <- logAR + d*(lgamma(alpha + beta + nOld[i]) - lgamma(alpha + beta + nNew[i]) ) + sum(lgamma(alpha + sNew[i,]) + lgamma(beta + nNew[i] - sNew[i,]) - lgamma(alpha + sOld[i,]) - lgamma(beta + nOld[i] - sOld[i,]))
			}
			if( log(runif(1)) < logAR ){
				reallocationAcceptanceRatio <- reallocationAcceptanceRatio + 1 
				z <- propZ
				s <- rep(0,K)
				sx <- array(data = 0, dim = c(K,d))
				for(k in 1:K){
					ind <- which(z == k)
					tmpV <- numeric(n)
					tmpM <- array(data = 0,dim = c(n,d))
					s[k] <- length(ind)
					if(s[k] > 0){
						tmpV[ind] <- rep(1,s[k])
						tmpM[ind,] <- x[ind,]
						sx[k,] <- colSums(array(x[ind,],dim = c(s[k],d)))
					}
					nMatrix[,k] <- rep(s[k],n) - tmpV
					sMatrix[,k,] <- matrix(sx[k,],nrow = n,ncol = d,byrow = TRUE) - tmpM
				}
			}
		}
		}
#		end of reallocation proposal 1
		if(mNumber == 'M2'){

#		reallocation proposal2
		if(K > 1){
			myPair <- sample(K,2,replace = FALSE)
			set1 <- which(z == myPair[1])
			set2 <- which(z == myPair[2])
			nOld <- c(length(set1),length(set2))
			propZ <- z
			if(nOld[1] > 0){
				randomSize <- 1 + floor(nOld[1]*runif(1))
				randomIndex <- set1[sample(nOld[1],randomSize,replace = FALSE)]
				propZ[randomIndex] <- rep(myPair[2],randomSize)
				newSet1 <- which(propZ == myPair[1])
				newSet2 <- which(propZ == myPair[2])
				nNew <- c(length(newSet1),length(newSet2))
				if(nOld[1] == 0){sOld1 <- rep(0,d)}else{sOld1 <- colSums(array(x[set1,],dim = c(nOld[1],d)))}
				if(nOld[2] == 0){sOld2 <- rep(0,d)}else{sOld2 <- colSums(array(x[set2,],dim = c(nOld[2],d)))}
				sOld <- rbind(sOld1,sOld2)
				if(nNew[1] == 0){sNew1 <- rep(0,d)}else{sNew1 <- colSums(array(x[newSet1,],dim = c(nNew[1],d)))}
				if(nNew[2] == 0){sNew2 <- rep(0,d)}else{sNew2 <- colSums(array(x[newSet2,],dim = c(nNew[2],d)))}
				sNew <- rbind(sNew1,sNew2)
				logAR <- log(nOld[1]) - log(nOld[2] + randomSize) -( lgamma(nOld[1] - randomSize + 1) + lgamma(nOld[2]+ randomSize + 1) - lgamma(nOld[2] + 1) - lgamma(nOld[1] + 1))
				for(i in 1:2){
					logAR <- logAR + heat*d*(lgamma(alpha + beta + nOld[i]) - lgamma(alpha + beta + nNew[i]) ) + heat*sum(lgamma(alpha + sNew[i,]) + lgamma(beta + nNew[i] - sNew[i,]) - lgamma(alpha + sOld[i,]) - lgamma(beta + nOld[i] - sOld[i,]))
				}
				logAR <- logAR + heat*sum(lgamma(gamma[myPair] + nNew)) - heat*sum(lgamma(gamma[myPair] + nOld)) 
				if( log(runif(1)) < logAR ){
					reallocationAcceptanceRatio2 <- reallocationAcceptanceRatio2 + 1 
					z <- propZ
					s <- rep(0,K)
					sx <- array(data = 0, dim = c(K,d))
					for(k in 1:K){
						ind <- which(z == k)
						tmpV <- numeric(n)
						tmpM <- array(data = 0,dim = c(n,d))
						s[k] <- length(ind)
						if(s[k] > 0){
							tmpV[ind] <- rep(1,s[k])
							tmpM[ind,] <- x[ind,]
							sx[k,] <- colSums(array(x[ind,],dim = c(s[k],d)))
						}
						nMatrix[,k] <- rep(s[k],n) - tmpV
						sMatrix[,k,] <- matrix(sx[k,],nrow = n,ncol = d,byrow = TRUE) - tmpM
					}
				}
			}
		}
		}
#		end of reallocation proposal 2
#		} #remove this
#		update the number of clusters
		if(mNumber == 'M3'){
		if ( runif(1) < birthProbs[K] ){
		# BIRTH MOVE
			birth = TRUE
			myPair <- c(sample(K,1),K + 1)
			# myPair[1]: ejecting component
			# myPair[2]: ejected component
			set1 <- which(z == myPair[1])
			nOld <- c(length(set1),0)
			ejectionAlpha <- toSolve(seq(0.0001,2,length = 500),nOld[1],pZERO);constLGAMMA <- 2*lgamma(ejectionAlpha) - lgamma(2*ejectionAlpha)

			propZ <- z
			myP <- rbeta(1,shape1 = ejectionAlpha,shape2 = ejectionAlpha)
			propZ[set1] <- myPair[sample(2,nOld[1], replace = TRUE,prob = c(1-myP,myP) )]
			newSet1 <- which(propZ == myPair[1])
			newSet2 <- which(propZ == myPair[2])
			nNew <- c(length(newSet1),length(newSet2))
			if(nOld[1] == 0){sOld1 <- rep(0,d)}else{sOld1 <- colSums(array(x[set1,],dim = c(nOld[1],d)))}
			if(nOld[2] == 0){sOld2 <- rep(0,d)}else{sOld2 <- colSums(array(x[set2,],dim = c(nOld[2],d)))}
			sOld <- rbind(sOld1,sOld2)
			if(nNew[1] == 0){sNew1 <- rep(0,d)}else{sNew1 <- colSums(array(x[newSet1,],dim = c(nNew[1],d)))}
			if(nNew[2] == 0){sNew2 <- rep(0,d)}else{sNew2 <- colSums(array(x[newSet2,],dim = c(nNew[2],d)))}
			sNew <- rbind(sNew1,sNew2)
			logAR <- log( (1 - birthProbs[K+1])/birthProbs[K] ) + constLGAMMA + lgamma(2*ejectionAlpha + nOld[1]) - lgamma(ejectionAlpha + nNew[1]) - lgamma(ejectionAlpha + nNew[2]) + heat*priorK[K+1] - heat*priorK[K]
			logAR <- logAR + heat*d*( lgamma(alpha + beta + nOld[1]) - lgamma(alpha + beta + nNew[1]) - lgamma(alpha + beta + nNew[2]) )
			logAR <- logAR + heat*sum( 
						  lgamma(alpha + sNew[1,]) + lgamma(beta + nNew[1] - sNew[1,]) 
						+ lgamma(alpha + sNew[2,]) + lgamma(beta + nNew[2] - sNew[2,]) 
						- lgamma(alpha + sOld[1,]) - lgamma(beta + nOld[1] - sOld[1,]) 
					)

#			for(i in 1:2){
#				logAR <- logAR + d*(lgamma(alpha + beta + nOld[i]) - lgamma(alpha + beta + nNew[i]) ) + sum(lgamma(alpha + sNew[i,]) + lgamma(beta + nNew[i] - sNew[i,]) - lgamma(alpha + sOld[i,]) - lgamma(beta + nOld[i] - sOld[i,]))
#			}
			logAR <- logAR + heat*sum(lgamma(gamma[myPair] + nNew)) - heat*sum(lgamma(gamma[myPair[1]] + nOld[1])) + heat*lgamma(sum(gamma[1:(K+1)])) - heat*lgamma(sum(gamma[1:K])) - heat*lgamma(sum(gamma[1:(K+1)]) + n) + heat*lgamma(sum(gamma[1:K]) + n) - heat*lgamma(gamma[K+1]) - heat*d*lbeta(alpha, beta)
		}else{
		# DEATH MOVE
			birth = FALSE
			myPair <- c(sample(K - 1,1),K)
			set1 <- which(z == myPair[1])
			set2 <- which(z == myPair[2])
			nOld <- c(length(set1),length(set2))
			propZ <- z
			propZ[c(set1,set2)] <- rep(myPair[1],sum(nOld))
			newSet1 <- which(propZ == myPair[1])
			nNew <- c(length(newSet1),0)
			ejectionAlpha <- toSolve(seq(0.0001,2,length = 500),nNew[1],pZERO);constLGAMMA <- 2*lgamma(ejectionAlpha) - lgamma(2*ejectionAlpha)
			if(nOld[1] == 0){sOld1 <- rep(0,d)}else{sOld1 <- colSums(array(x[set1,],dim = c(nOld[1],d)))}
			if(nOld[2] == 0){sOld2 <- rep(0,d)}else{sOld2 <- colSums(array(x[set2,],dim = c(nOld[2],d)))}
			sOld <- rbind(sOld1,sOld2)
			if(nNew[1] == 0){sNew1 <- rep(0,d)}else{sNew1 <- colSums(array(x[newSet1,],dim = c(nNew[1],d)))}
			sNew2 <- rep(0,d)
			sNew <- rbind(sNew1,sNew2)
			logAR <- - (log( (1 - birthProbs[K])/birthProbs[K - 1] ) + constLGAMMA + lgamma(2*ejectionAlpha + nNew[1]) - lgamma(ejectionAlpha + nOld[1]) - lgamma(ejectionAlpha + nOld[2]))+ heat*priorK[K - 1] - heat*priorK[K]


			logAR <- logAR + heat*d*( lgamma(alpha + beta + nOld[1]) + lgamma(alpha + beta + nOld[2]) - lgamma(alpha + beta + nNew[1]) )
			logAR <- logAR + heat*sum( 
						  lgamma(alpha + sNew[1,]) + lgamma(beta + nNew[1] - sNew[1,]) 
						- lgamma(alpha + sOld[1,]) - lgamma(beta + nOld[1] - sOld[1,]) 
						- lgamma(alpha + sOld[2,]) - lgamma(beta + nOld[2] - sOld[2,])
					 )


#			for(i in 1:2){
#				logAR <- logAR + d*(lgamma(alpha + beta + nOld[i]) - lgamma(alpha + beta + nNew[i]) ) + sum(lgamma(alpha + sNew[i,]) + lgamma(beta + nNew[i] - sNew[i,]) - lgamma(alpha + sOld[i,]) - lgamma(beta + nOld[i] - sOld[i,]))
#			}
			logAR <- logAR + heat*lgamma(gamma[myPair[1]] + nNew[1]) - heat*sum(lgamma(gamma[myPair] + nOld)) + heat*lgamma(sum(gamma[1:(K-1)])) - heat*lgamma(sum(gamma[1:K])) - heat*lgamma(sum(gamma[1:(K-1)]) + n) + heat*lgamma(sum(gamma[1:K]) + n) + heat*lgamma(gamma[K]) + heat*d*lbeta(alpha, beta)
			
		}

		if( log(runif(1)) < logAR ){
			reallocationAcceptanceRatio3 <- reallocationAcceptanceRatio3 + 1 
			z <- propZ
			if(birth == TRUE){K <- K + 1}else{K <- K - 1}
			s <- rep(0,K)
			sx <- array(data = 0, dim = c(K,d))
			for(k in 1:K){
				ind <- which(z == k)
				tmpV <- numeric(n)
				tmpM <- array(data = 0,dim = c(n,d))
				s[k] <- length(ind)
				if(s[k] > 0){
					tmpV[ind] <- rep(1,s[k])
					tmpM[ind,] <- x[ind,]
					sx[k,] <- colSums(array(x[ind,],dim = c(s[k],d)))
				}
				nMatrix[,k] <- rep(s[k],n) - tmpV
				sMatrix[,k,] <- matrix(sx[k,],nrow = n,ncol = d,byrow = TRUE) - tmpM
			}
		}
		}
#		end of update for the number of clusters


		if(mNumber == 'M4'){

#	reallocation proposal from conditional probs

			if(K > 1){
			myPair <- sample(K,2,replace = FALSE)
			propZ <- z
			set1 <- which(z == myPair[1])
			set2 <- which(z == myPair[2])
			union12 <- c(set1,set2)
			processed <- c()
			notProcessed <- setdiff(union12,processed)
			n1New <- n2New <- 0
			s1New <- s2New <- numeric(d)
			n1 <- n2 <- 0
			s1 <- s2 <- numeric(d)
			j1 <- myPair[1]
			j2 <- myPair[2]
			alreadyInJ1 <- alreadyInJ2 <- c()
			proposalRatio <- 0
			checkCondition <- TRUE
			while ( length(notProcessed) > 0 ){
				i <- notProcessed[1]
				processed <- c(processed,i)
				# propose reallocation
				u <- log(gamma[j1] + n1New) - log(gamma[j2] + n2New) + sum( lgamma(alpha + x[i,] + s1New) + lgamma(beta + 1 + n1New - s1New - x[i,]) ) - d*sum(lgamma(alpha+beta+n1New+1)) + sum( lgamma(alpha + s2New) + lgamma(beta + n2New - s2New) ) - d*sum(lgamma(alpha+beta+n2New)) - sum( lgamma(alpha + s1New) + lgamma(beta + n1New - s1New) ) + d*sum(lgamma(alpha+beta+n1New)) - sum( lgamma(alpha + x[i,] + s2New) + lgamma(beta + 1 + n2New - s2New - x[i,]) ) + d*sum(lgamma(alpha+beta+n2New+1))
				u <- exp(u)
				proposalProbabities <- u/(1 + u)
				if( is.finite(proposalProbabities) == FALSE ){checkCondition <- FALSE; proposalProbabities <- 0.5}
				if( runif(1) < proposalProbabities ){ 
					propZ[i] <- j1
					proposalRatio <- proposalRatio - log(proposalProbabities) 
					n1New <- n1New + 1
					s1New <- s1New + x[i,]
				}else{ 
					propZ[i] <- j2
					proposalRatio <- proposalRatio - log(1 - proposalProbabities) 
					n2New <- n2New + 1
					s2New <- s2New + x[i,]
				}
				# compute inverse move
				u <- log(gamma[j1] + n1) - log(gamma[j2] + n2) + sum( lgamma(alpha + x[i,] + s1) + lgamma(beta + 1 + n1 - s1 - x[i,]) ) - d*sum(lgamma(alpha+beta+n1+1)) + sum( lgamma(alpha + s2) + lgamma(beta + n2 - s2) ) - d*sum(lgamma(alpha+beta+n2)) - sum( lgamma(alpha + s1) + lgamma(beta + n1 - s1) ) + d*sum(lgamma(alpha+beta+n1)) - sum( lgamma(alpha + x[i,] + s2) + lgamma(beta + 1 + n2 - s2 - x[i,]) ) + d*sum(lgamma(alpha+beta+n2+1))
				u <- exp(u)
				u <- u/(1 + u)
				if(z[i] == j1){ 
					proposalRatio <- proposalRatio + log(u) 
					n1 <- n1 + 1
					s1 <- s1 + x[i,]
				}else{ 
					proposalRatio <- proposalRatio + log(1 - u) 
					n2 <- n2 + 1
					s2 <- s2 + x[i,]
				}
				
				notProcessed <- setdiff(union12,processed)
			}

			newSet1 <- which(propZ == myPair[1])
			newSet2 <- which(propZ == myPair[2])
			nNew <- c(length(newSet1),length(newSet2))
			nOld <- c(length(set1),length(set2))
			if(nOld[1] == 0){sOld1 <- rep(0,d)}else{sOld1 <- colSums(array(x[set1,],dim = c(nOld[1],d)))}
			if(nOld[2] == 0){sOld2 <- rep(0,d)}else{sOld2 <- colSums(array(x[set2,],dim = c(nOld[2],d)))}
			sOld <- rbind(sOld1,sOld2)
			if(nNew[1] == 0){sNew1 <- rep(0,d)}else{sNew1 <- colSums(array(x[newSet1,],dim = c(nNew[1],d)))}
			if(nNew[2] == 0){sNew2 <- rep(0,d)}else{sNew2 <- colSums(array(x[newSet2,],dim = c(nNew[2],d)))}
			sNew <- rbind(sNew1,sNew2)
			logAR <- proposalRatio
			for (i in 1:2){
				logAR <- logAR + d*(lgamma(alpha + beta + nOld[i]) - lgamma(alpha + beta + nNew[i]) ) + sum(lgamma(alpha + sNew[i,]) + lgamma(beta + nNew[i] - sNew[i,]) - lgamma(alpha + sOld[i,]) - lgamma(beta + nOld[i] - sOld[i,]))
			}
			logAR <- logAR + sum(lgamma(gamma[myPair] + nNew)) - sum(lgamma(gamma[myPair] + nOld)) 
			if( (log(runif(1)) < logAR) && (checkCondition == TRUE) ){
				reallocationAcceptanceRatio4 <- reallocationAcceptanceRatio4 + 1 
				z <- propZ
				s <- rep(0,K)
				sx <- array(data = 0, dim = c(K,d))
				for(k in 1:K){
					ind <- which(z == k)
					tmpV <- numeric(n)
					tmpM <- array(data = 0,dim = c(n,d))
					s[k] <- length(ind)
					if(s[k] > 0){
						tmpV[ind] <- rep(1,s[k])
						tmpM[ind,] <- x[ind,]
						sx[k,] <- colSums(array(x[ind,],dim = c(s[k],d)))
					}
					nMatrix[,k] <- rep(s[k],n) - tmpV
					sMatrix[,k,] <- matrix(sx[k,],nrow = n,ncol = d,byrow = TRUE) - tmpM
				}
			}

			}
#	end of reallocation proposal from conditional probs

			}







		if(iter %% thinning == 0){
			kValues <- c(kValues,K)
			if(iter > burn){
				cat(K,"\n",file=conK)
				cat(z,"\n",file=conZ)
				if (heat == 1){ #produce parameter values only for cold chain
					for(k in 1:K){
						myIndex <- which(z == k)
						s[k] <- length(myIndex)
						if(s[k] > 1){
							sx[k,] <- colSums(x[myIndex,])
						}else{
							if(s[k] == 1){
								sx[k,] <- x[myIndex,]
							}else{
								sx[k,] <- rep(0,d)
							}

						}
						theta[k,] <- rbeta(d, shape1 = alpha + sx[k,], shape2 = beta + s[k] - sx[k,])
					}
					p <- myDirichlet(gamma[1:K] + s[1:K])
					cat(theta,"\n",file=conTheta)
					cat(p,"\n",file=conP)
				}		
			}
		}
		if(iter %% (m/100) == 0){
			#if(length(kValues) > 0){
			#	plot(kValues)
			#}
#			if(length(kValues) > 0 ){points(rep(iter,K),p,col = 1:K,pch = K)}
			cat(paste("                             ",100*round(iter/m,3),"% completed. ", 100*round(lmm*reallocationAcceptanceRatio/iter,3),"%, ",100*round(lmm*reallocationAcceptanceRatio2/iter,3),"%, ",100*round(lmm*reallocationAcceptanceRatio4/iter,3),"%, ",100*round(lmm*reallocationAcceptanceRatio3/iter,3),"%",sep=""),"\n");}

		# label switching move
#		perm <- sample(K,K,replace = FALSE)
#		if(K < Kmax){
#			perm <- c(perm,(K+1):Kmax)
#		}
#		z <- perm[z]
#		nMatrix <- nMatrix[,perm]
#		sMatrix <- sMatrix[,perm,] 

	}
	close(conK)
	close(conTheta)
	close(conP)
	close(conZ)
	cat(paste("Allocation Sampler finished."),"\n")
	cat("\n")
	#if(K < 1){
		if(LS == TRUE){
		kFile <- read.table("K.txt")[,1]
		cat("\n")
		print(table(kFile)/length(kFile))
		K <- as.numeric(names(table(kFile))[order(table(kFile),decreasing=TRUE)[1]])
		cat("\n")
		cat(paste("Most probable model: K = ",K," with P(K = ",K,") = ",max(table(kFile)/length(kFile)),sep=""),"\n")		
		cat("\n")
	if(missing(reorderModels)==TRUE){reorderModels <- 'onlyMAP'}
	if(reorderModels == 'onlyMAP'){kRange <- K}else{
		kRange <- as.numeric(names(table(kFile)/length(kFile)))
	}
	for ( K in kRange ){
		if (K > 1){
			cat(paste("Dealing with Label Switching for K =",K),"\n")
			tt <- read.table("theta.varK.txt")
			index <- which(kFile == K)
			d <- dim(x)[2]
			J <- d + 1
			m <- length(index)
			mcmc <- array(data = NA, dim = c(m,K,J))
			Kmax <- (dim(tt)/d)[2]
			for (j in 1:d){
				for(k in 1:K){
					mcmc[,k,j] <- tt[index,(j-1)*Kmax + k]
				}
			}
			conP = file(p.file,open = "r")
			i <- 0
			j <- 0
			while (length(oneLine <- readLines(conP, n = 1, warn = FALSE)) > 0) {
				i <- i + 1
				if(kFile[i] == K){
					j <- j + 1
					mcmc[j,,J] <- as.numeric(strsplit(oneLine,split = " ")[[1]])
				}
			}
			close(conP)
			allocations <- as.matrix(read.table("z.varK.txt",as.is = TRUE)[index,])
			iter <- 1
			ll <- complete.loglikelihood(x,allocations[iter,],mcmc[iter,,])
			maxLL <- ll
			maxIter <- iter
			for (iter in 1:m){
				ll <- complete.loglikelihood(x,allocations[iter,],mcmc[iter,,])
				if(ll > maxLL){maxLL <- ll;maxIter <- iter;cat(paste("Found new Complete MLE: ", ll,sep=""),"\n")}
			}
			# classification probs
			pMatrix <- array(data = NA, dim = c(m,n,K))
			l <- numeric(K)
			for (iter in 1:m){
				p <- mcmc[iter,,J]
				theta <- array( mcmc[iter,,1:(J-1)],dim = c(K,J-1) )
				for(i in 1:n){
					for(k in 1:K){
						l[k] <- log(p[k]) + sum(x[i,]*log(theta[k,]) + (1-x[i,])*log(1-theta[k,]))
					}
					for(k in 1:K){
						pMatrix[iter,i,k] <- 1/sum(exp(l - l[k]))
						if(is.na(pMatrix[iter,i,k]) == TRUE){pMatrix[iter,i,k] = 0}
					}
				}
				if(iter %% 1000 == 0){cat(paste(" classification probs: ",100*round(iter/m,3),"% completed",sep=""),"\n");}
			}


			if(missing(z.true)==TRUE){
				ls <- label.switching( method = c("STEPHENS","ECR","ECR-ITERATIVE-1"),
							zpivot = allocations[maxIter,], z = allocations,K = K, complete = complete.loglikelihood, data = x,
							prapivot = mcmc[maxIter,,], mcmc = mcmc, p = pMatrix)
			}else{
				ls <- label.switching( method = c("STEPHENS","ECR","ECR-ITERATIVE-1"),
							zpivot = allocations[maxIter,], z = allocations,K = K, complete = complete.loglikelihood, data = x,
							prapivot = mcmc[maxIter,,], mcmc = mcmc, p = pMatrix,groundTruth = z.true)
			}
			reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"ECR")$output
			write.table(reordered.mcmc, file = paste("reorderedMCMC-ECR.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"ECR-ITERATIVE-1")$output
			write.table(reordered.mcmc, file = paste("reorderedMCMC-ECR-ITERATIVE1.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"STEPHENS")$output
			write.table(reordered.mcmc, file = paste("reorderedMCMC-STEPHENS.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			write.table(mcmc, file = paste("rawMCMC.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			write.table(file = paste("reorderedSingleBestClusterings.mapK.",K,".txt",sep=""),t(ls$clusters[c(1,2,3),]),row.names = paste("z",1:n,sep="."))

			#reordering allocations
			allocationsECR <- allocationsKL <- allocationsECR.ITERATIVE1 <- allocations
			for (i in 1:m){
				myPerm <- order(ls$permutations$"ECR"[i,])
				allocationsECR[i,] <- myPerm[allocations[i,]]
				myPerm <- order(ls$permutations$"STEPHENS"[i,])
				allocationsKL[i,] <- myPerm[allocations[i,]]
				myPerm <- order(ls$permutations$"ECR-ITERATIVE-1"[i,])
				allocationsECR.ITERATIVE1[i,] <- myPerm[allocations[i,]]

			}
			write.table(allocationsECR, file = paste("z.ECR.mapK.",K,".txt",sep=""))
			write.table(allocationsKL, file = paste("z.KL.mapK.",K,".txt",sep=""))
			write.table(allocationsECR.ITERATIVE1, file = paste("z.ECR-ITERATIVE1.mapK.",K,".txt",sep=""))
			#file.remove("p.collapsed.txt")
			#file.remove("theta.collapsed.txt")
			cat(paste0("raw MCMC parameters for most probable K written to: \'rawMCMC.mapK.",K,".txt\' "),"\n")
			cat(paste0("raw MCMC latent allocations for most probable K written to: \'z.mapK.",K,".txt\' "),"\n")
			cat(paste0("reordered MCMC output written to: "),"\n")
			cat(paste0("     (Method 1):     \'reorderedMCMC-ECR.mapK.",K,".txt\'"),"\n")
			cat(paste0("     (Method 2):     \'reorderedMCMC-ECR-ITERATIVE1.mapK.",K,".txt\'"),"\n")
			cat(paste0("     (Method 3):     \'reorderedMCMC-STEPHENS.mapK.",K,".txt\'"),"\n")
			cat(paste0("reordered single best clusterings written to: \'reorderedSingleBestClusterings.mapK.",K,".txt\' "),"\n")
			cat(paste0("reordered MCMC latent allocations for most probable K written to: "),"\n")
			cat(paste0("     (Method 1):     \'z.ECR.mapK.",K,".txt\'"),"\n")
			cat(paste0("     (Method 2):     \'z.KL.mapK.",K,".txt\'"),"\n")
			cat(paste0("     (Method 3):     \'z.ECR-ITERATIVE1.mapK.",K,".txt\'"),"\n")

		}

	}
	}
	setwd("../")
	cat(paste0("back to working directory now."),"\n")
}

coupledMetropolis <- function(Kmax, nChains,heats,binaryData,outPrefix,ClusterPrior,m, alpha, beta, gamma, z.true, ejectionAlpha){
	if(missing(nChains) == TRUE){stop(cat(paste("    [ERROR]: number of chains not provided."), "\n"))}
	if(missing(heats) == TRUE){
		heats <- seq(1,0.1,length = nChains)
		}else{
		if(heats[1] != 1){stop(cat(paste("    [ERROR]: `heats[1]` should be equal to one."), "\n"))}
	}
	if(length(heats) != nChains){
		stop(cat(paste("    [ERROR]: `length(heats)` should be equal to `nChains`."), "\n"))
	}
	if(missing(ClusterPrior) == TRUE){ClusterPrior <- 'poisson'}
	if (missing(alpha)) {alpha <- 1}
	if (missing(beta)) {beta <- 1}
	if (missing(binaryData)) {stop(cat(paste("    [ERROR]: data is missing."), "\n"))}
	if (missing(gamma)) {gamma <- rep(1,Kmax)}
	if (length(table(gamma)) > 1){
		stop(cat(paste("    [ERROR]: Dirichlet prior parameters should be the same."), "\n"))
	}
	d <- dim(binaryData)[2]; n <- dim(binaryData)[1]; priorK <- numeric(Kmax)
	if(ClusterPrior == "uniform"){
		priorK <- rep(log(1/Kmax),Kmax)
	}
	if(ClusterPrior == "poisson"){
		denom <- log(ppois(Kmax, lambda = 1, lower.tail = TRUE) - dpois(0,lambda = 1))
		for (k in 1:Kmax){
			priorK[k] <- dpois(k,lambda = 1, log = TRUE) - denom 
		}
	}
	if(missing(ejectionAlpha) == TRUE){ejectionAlpha <- 0.2}
	x <- binaryData
	currentZ <- array(data = NA, dim = c(nChains,n))
	currentK <- numeric(nChains)
	# check if dir exists and stop in this case
	myCheck <- match(outPrefix,list.files())
	if(is.na(myCheck) == FALSE){
		stop(cat(paste("    [ERROR]: directory exists, please provide different name to outPrefix."), "\n"))
	}
	if( missing(z.true) ){
		z.true <- NULL
	}



	if(nChains < 2){
		cat(paste("            Only 1 chain? Well."),"\n")
		dir.create(outPrefix)
		myHeat <- 1
		allocationSamplerBinMix( alpha = alpha, beta = beta, gamma = gamma, m = 10*m, burn= 9, data = binaryData, 
					thinning = 10,Kmax = Kmax, ClusterPrior = ClusterPrior,ejectionAlpha = ejectionAlpha, 
					outputDir = outPrefix,Kstart = 1,heat=myHeat,metropolisMoves = c('M1','M2','M3','M4'),LS = FALSE)}
	else{

		registerDoParallel(cores = nChains)
		#registerDoMC(nChains)
		outputDirs <- paste0(outPrefix,1:nChains)
		temperatures <- heats
		myChain <- 1
		for(i in 1:nChains){
			cat(paste0("            Create temporary directory: \'",outPrefix,i,"\'."), "\n")
		}
		cat(paste0("    [NOTE]: screen output from multiple threads is redirected to \'",outPrefix,"-sampler.log\'."), "\n")
		sink(paste0(outPrefix,'-sampler.log'))
		foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
			outDir <- outputDirs[myChain]
			dir.create(outDir)
			myHeat <- temperatures[myChain]
			allocationSamplerBinMix( alpha = alpha, beta = beta, gamma = gamma, m = 10, burn= 9, data = binaryData, 
						thinning = 1,Kmax = Kmax, ClusterPrior = ClusterPrior,ejectionAlpha = ejectionAlpha, 
						outputDir = outDir,Kstart=1,heat=myHeat,metropolisMoves='M3',LS = FALSE)
		}



		ITERATIONS <- m
		sampledK <- array(data = 0, dim = c(ITERATIONS,nChains))
		dir.create(outPrefix)
		k.file <- paste0(outPrefix,"/K.txt")
		theta.file <- paste0(outPrefix,"/theta.varK.txt")
		z.file <- paste0(outPrefix,"/z.varK.txt")
		p.file <- paste0(outPrefix,"/p.varK.txt")
		conK = file(k.file,open = "w")
		conTheta = file(theta.file,open = "w")
		conZ = file(z.file,open = "w")
		conP = file(p.file,open = "w")
		ar <- 0
		metMoves <- vector('list',length = nChains)
		metMoves[[1]] <- c('M1','M2','M3','M4')
		for(j in 2:nChains){metMoves[[j]] <- c('M2','M3')}
		localAR <- 0	#every 10 iterations
		for(iter in 1:ITERATIONS){
			for(myChain in 1:nChains){
				currentZ[myChain,] <- as.numeric(read.table(paste0(outPrefix,myChain,'/z.varK.txt'))[1,])
				currentK[myChain] <- read.table(paste0(outPrefix,myChain,'/K.txt'))[1,]
			}
			sampledK[iter,] <- currentK
			myPair <- sample(nChains,2,replace = FALSE)
			j1 <- myPair[1]
			j2 <- myPair[2]
			z1 <- currentZ[j1,]
			z2 <- currentZ[j2,]
			K1 <- currentK[j1]
			K2 <- currentK[j2]

			s1 <- rep(0,K1)
			sx1 <- array(data = 0, dim = c(K1,d))
			for(k in 1:K1){
				ind <- which(z1 == k)
				s1[k] <- length(ind)
				if(s1[k] > 0){
					sx1[k,] <- colSums(array(x[ind,],dim = c(s1[k],d)))
				}
			}
			s2 <- rep(0,K2)
			sx2 <- array(data = 0, dim = c(K2,d))
			for(k in 1:K2){
				ind <- which(z2 == k)
				s2[k] <- length(ind)
				if(s2[k] > 0){
					sx2[k,] <- colSums(array(x[ind,],dim = c(s2[k],d)))
				}
			}
			#compute log-posterior for 1st
			log.posterior <- priorK[K1] + sum(lgamma( gamma[1:K1] + s1 )) - d*sum(lgamma(alpha+beta+s1)) - lgamma(n + sum(gamma[1:K1]))
			for(k in 1:K1){
				log.posterior <- log.posterior + sum(lgamma(alpha + sx1[k,]) + lgamma(beta + s1[k] - sx1[k,]))
			}
			#prior constants1:
			log.posterior <- log.posterior + lgamma(sum(gamma[1:K1])) - sum(lgamma(gamma[1:K1])) - K1*d*lbeta(alpha, beta)
			logAR <- (temperatures[j2] - temperatures[j1])*log.posterior
			#compute log-posterior for 2nd
			log.posterior <- priorK[K2] + sum(lgamma( gamma[1:K2] + s2 )) - d*sum(lgamma(alpha+beta+s2)) - lgamma(n + sum(gamma[1:K2]))
			for(k in 1:K2){
				log.posterior <- log.posterior + sum(lgamma(alpha + sx2[k,]) + lgamma(beta + s2[k] - sx2[k,]))
			}
			#prior constants2:
			log.posterior <- log.posterior + lgamma(sum(gamma[1:K2])) - sum(lgamma(gamma[1:K2])) - K2*d*lbeta(alpha, beta)
			logAR <- logAR + (temperatures[j1] - temperatures[j2])*log.posterior
			#cat(paste('logAR =',logAR),'\n')
			if( log(runif(1)) < logAR ){
				ar <- ar + 1
				localAR <- localAR + 1
				currentZ[j1,] <- z2
				currentZ[j2,] <- z1
				currentK[j1] <- K2
				currentK[j2] <- K1
				#cat(paste('switching.'),'\n')
			}
			myChain <- 1
			foreach(myChain=1:nChains, .export=ls(envir=globalenv()) ) %dopar% {
				outDir <- outputDirs[myChain]
				myHeat <- temperatures[myChain]
				Kstart <- currentK[myChain]
				zStart <- currentZ[myChain,]
				allocationSamplerBinMix( alpha = alpha, beta = beta, gamma = gamma, m = 10, burn= 9, data = binaryData, 
							thinning = 1,Kmax = Kmax, ClusterPrior = ClusterPrior,ejectionAlpha = ejectionAlpha, 
							outputDir = outDir,Kstart=Kstart,zStart = zStart, heat=myHeat,metropolisMoves =  metMoves[[myChain]],LS = FALSE)
			}

			if(iter %% (m/100) == 0){
				matplot(sampledK[1:iter,],type = "l",lty = 1,lwd = 2, col = topo.colors(nChains))
				legend('topleft',paste0('f(z,K|data)^{',round(heats,3),'}'),lty = 1, lwd = 2, col = topo.colors(nChains))
				write(paste0(100*iter/m,'% completed. Chain switching acceptance rate: ',100*round(ar/iter,3),'%.'),stderr())
			}
			kk <- as.numeric(read.table(paste0(outPrefix,"1/K.txt"))[1,])
			cat(kk,"\n",file=conK)
			theta <- as.numeric(read.table(paste0(outPrefix,"1/theta.varK.txt"))[1,])
			cat(theta,"\n",file=conTheta)
			z <- as.numeric(read.table(paste0(outPrefix,"1/z.varK.txt"))[1,])
			cat(z,"\n",file=conZ)
			p <- as.numeric(read.table(paste0(outPrefix,"1/p.varK.txt"))[1,])
			cat(p,"\n",file=conP)		

		}
		close(conK)
		close(conTheta)
		close(conP)
		close(conZ)
		write.table(sampledK, file = paste0(outPrefix,"/K.allChains.txt"), col.names = paste0('chain.',1:nChains),quote=FALSE,row.names = FALSE)
		for(k in 1:nChains){
			unlink(paste0(outPrefix,k), recursive=TRUE)
		}
		sink()
		stopImplicitCluster()
	}
	if(is.null(z.true) == TRUE){
		dealWithLabelSwitching(outDir = outPrefix, binaryData = binaryData)
	}else{
		dealWithLabelSwitching(outDir = outPrefix, binaryData = binaryData, z.true = z.true)
	}
}



dealWithLabelSwitching <- function(outDir,reorderModels, binaryData,z.true){
	setwd(outDir)
	x <- binaryData
	kFile <- read.table("K.txt")[,1]
	theta.file <- "theta.varK.txt"
	z.file <- "z.varK.txt"
	p.file <- "p.varK.txt"
	cat("\n")
	print(table(kFile)/length(kFile))
	K <- as.numeric(names(table(kFile))[order(table(kFile),decreasing=TRUE)[1]])
	cat("\n")
	cat(paste("Most probable model: K = ",K," with P(K = ",K,") = ",max(table(kFile)/length(kFile)),sep=""),"\n")		
	cat("\n")
	if(missing(reorderModels)==TRUE){reorderModels <- 'onlyMAP'}
	if(reorderModels == 'onlyMAP'){kRange <- K}else{
		kRange <- as.numeric(names(table(kFile)/length(kFile)))
	}
	for ( K in kRange ){

		tt <- read.table("theta.varK.txt")
		index <- which(kFile == K)
		d <- dim(x)[2]
		n <- dim(x)[1]
		J <- d + 1
		m <- length(index)
		mcmc <- array(data = NA, dim = c(m,K,J))
		Kmax <- (dim(tt)/d)[2]
		for (j in 1:d){
			for(k in 1:K){
				mcmc[,k,j] <- tt[index,(j-1)*Kmax + k]
			}
		}
		conP = file(p.file,open = "r")
		i <- 0
		j <- 0
		while (length(oneLine <- readLines(conP, n = 1, warn = FALSE)) > 0) {
			i <- i + 1
			if(kFile[i] == K){
				j <- j + 1
				mcmc[j,,J] <- as.numeric(strsplit(oneLine,split = " ")[[1]])
			}
		}
		close(conP)
		if (K > 1){
			cat(paste("Dealing with Label Switching for K =",K),"\n")
			allocations <- as.matrix(read.table("z.varK.txt",as.is = TRUE)[index,])
			iter <- 1
			ll <- complete.loglikelihood(x,allocations[iter,],mcmc[iter,,])
			maxLL <- ll
			maxIter <- iter
			for (iter in 1:m){
				ll <- complete.loglikelihood(x,allocations[iter,],mcmc[iter,,])
				if(ll > maxLL){maxLL <- ll;maxIter <- iter;cat(paste("Found new Complete MLE: ", ll,sep=""),"\n")}
			}
			# classification probs
			#####OLD 
#			ptm <- proc.time()
#			pMatrix <- array(data = NA, dim = c(m,n,K))
#			l <- numeric(K)
#			for (iter in 1:m){
#				p <- mcmc[iter,,J]
#				theta <- array( mcmc[iter,,1:(J-1)],dim = c(K,J-1) )
#				for(i in 1:n){
#					for(k in 1:K){
#						l[k] <- log(p[k]) + sum(x[i,]*log(theta[k,]) + (1-x[i,])*log(1-theta[k,]))
#					}
#					for(k in 1:K){
#						pMatrix[iter,i,k] <- 1/sum(exp(l - l[k]))
#						if(is.na(pMatrix[iter,i,k]) == TRUE){pMatrix[iter,i,k] = 0}
#					}
#				}
#				if(iter %% 1000 == 0){cat(paste(" classification probs: ",100*round(iter/m,3),"% completed",sep=""),"\n");}
#			}
#			cat(paste0("proc.time for classification probabalities1: ", round(as.numeric((proc.time() - ptm)[3]),2)),"\n")
			#####NEW 
			ptm <- proc.time()
			pMatrix <- array(data = NA, dim = c(m,n,K))
			l <- array(data = 0, dim = c(n,K))
			for (iter in 1:m){
				LOG.P <- log(mcmc[iter,,J])
				LOG.THETA <- array(log(mcmc[iter,,1:(J-1)]),dim = c(K,J-1) )
				LOG.1_MINUS_THETA <- array(log(1 - mcmc[iter,,1:(J-1)]),dim = c(K,J-1) )
					for(k in 1:K){
						l[,k] <- LOG.P[k] + rowSums( x * matrix(LOG.THETA[k,], nrow = n, ncol = J - 1, byrow = TRUE) + (1-x)*matrix(LOG.1_MINUS_THETA[k,], nrow = n, ncol = J - 1, byrow = TRUE))
					}
					for(k in 1:K){
						pMatrix[iter,,k] <- apply(l,1,function(y){return(1/sum(exp(y - y[k])))} )
						ind <- which(is.na(pMatrix[iter,,k]) == TRUE) 
						nInd <- length(ind)
						if(nInd > 0){pMatrix[iter,ind,k] <- rep(0,nInd)}
					}
				
				if(iter %% 1000 == 0){cat(paste(" classification probs: ",100*round(iter/m,3),"% completed",sep=""),"\n");}
			}

			cat(paste0("proc.time for classification probabilities: ", round(as.numeric((proc.time() - ptm)[3]),2)),"\n")
			if(missing(z.true)==TRUE){
				ls <- label.switching( method = c("STEPHENS","ECR","ECR-ITERATIVE-1"),
							zpivot = allocations[maxIter,], z = allocations,K = K, complete = complete.loglikelihood, data = x,
							prapivot = mcmc[maxIter,,], mcmc = mcmc, p = pMatrix)
			}else{
				ls <- label.switching( method = c("STEPHENS","ECR","ECR-ITERATIVE-1"),
							zpivot = allocations[maxIter,], z = allocations,K = K, complete = complete.loglikelihood, data = x,
							prapivot = mcmc[maxIter,,], mcmc = mcmc, p = pMatrix,groundTruth = z.true)
			}
			reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"ECR")$output
			write.table(reordered.mcmc, file = paste("reorderedMCMC-ECR.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"ECR-ITERATIVE-1")$output
			write.table(reordered.mcmc, file = paste("reorderedMCMC-ECR-ITERATIVE1.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			reordered.mcmc<-permute.mcmc(mcmc,ls$permutations$"STEPHENS")$output
			write.table(reordered.mcmc, file = paste("reorderedMCMC-STEPHENS.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			write.table(mcmc, file = paste("rawMCMC.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			write.table(file = paste("reorderedSingleBestClusterings.mapK.",K,".txt",sep=""),t(ls$clusters[c(1,2,3),]),row.names = paste("z",1:n,sep="."))

			#reordering allocations 
			allocationsECR <- allocationsKL <- allocationsECR.ITERATIVE1 <- allocations
			for (i in 1:m){
				myPerm <- order(ls$permutations$"ECR"[i,])
				allocationsECR[i,] <- myPerm[allocations[i,]]
				myPerm <- order(ls$permutations$"STEPHENS"[i,])
				allocationsKL[i,] <- myPerm[allocations[i,]]
				myPerm <- order(ls$permutations$"ECR-ITERATIVE-1"[i,])
				allocationsECR.ITERATIVE1[i,] <- myPerm[allocations[i,]]

			}
			MeanReorderedpMatrix <- array(data = 0, dim = c(n,K))    # define object that will contain the classification probs
			for (i in 1:m){
				myPerm <- ls$permutations$"ECR"[i,]   # this is the permutation of labels for iteration i according to ECR algorithm
				MeanReorderedpMatrix <- MeanReorderedpMatrix + pMatrix[i, ,myPerm]   # apply myPerm to the columns of pMatrix for given iteration and add the permuted matrix to MeanReorderedpMatrix
			}
			MeanReorderedpMatrix <- MeanReorderedpMatrix/m    # this is the final estimate of classification probabilities. 
			write.csv(MeanReorderedpMatrix, file = paste0("classificationProbabilities.mapK.",K,".csv"), row.names = FALSE)


			write.table(allocationsECR, file = paste("z.ECR.mapK.",K,".txt",sep=""))
			write.table(allocationsKL, file = paste("z.KL.mapK.",K,".txt",sep=""))
			write.table(allocationsECR.ITERATIVE1, file = paste("z.ECR-ITERATIVE1.mapK.",K,".txt",sep=""))
			#file.remove("p.collapsed.txt")
			#file.remove("theta.collapsed.txt")
			cat(paste0("raw MCMC parameters for most probable K written to: \'rawMCMC.mapK.",K,".txt\' "),"\n")
			cat(paste0("raw MCMC latent allocations for most probable K written to: \'z.mapK.",K,".txt\' "),"\n")
			cat(paste0("reordered MCMC output written to: "),"\n")
			cat(paste0("     (Method 1):     \'reorderedMCMC-ECR.mapK.",K,".txt\'"),"\n")
			cat(paste0("     (Method 2):     \'reorderedMCMC-ECR-ITERATIVE1.mapK.",K,".txt\'"),"\n")
			cat(paste0("     (Method 3):     \'reorderedMCMC-STEPHENS.mapK.",K,".txt\'"),"\n")
			cat(paste0("reordered single best clusterings written to: \'reorderedSingleBestClusterings.mapK.",K,".txt\' "),"\n")
			cat(paste0("reordered MCMC latent allocations for most probable K written to: "),"\n")
			cat(paste0("     (Method 1):     \'z.ECR.mapK.",K,".txt\'"),"\n")
			cat(paste0("     (Method 2):     \'z.KL.mapK.",K,".txt\'"),"\n")
			cat(paste0("     (Method 3):     \'z.ECR-ITERATIVE1.mapK.",K,".txt\'"),"\n")

		}else{
			cat(paste0('[NOTE]:    Most probable model corresponds to 1 cluster so the label-switching algorithms are not applied.',"\n"))
			write.table(mcmc, file = paste("MCMC.mapK.",K,".txt",sep=""),col.names = c(paste(rep(paste(expression(theta),1:K,sep="."),d),rep(1:d,each=K),sep="-"),paste('p',1:K,sep=".")))
			cat(paste0("MCMC output corresponding to most probable model (K = 1) written to: \'", paste("MCMC.mapK.",K,".txt",sep=""),"\'"),"\n")


		}

	}
	setwd("../")

}






