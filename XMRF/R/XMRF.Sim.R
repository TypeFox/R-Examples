XMRF.Sim <-
function(n=100, p=50, model="LPGM", graph.type="scale-free"){
	mydata <- list()
	
	# Generate a network first
	B <- simGraph(p=p, type=graph.type)
	mydata$B <- B
	
	if(model == "LPGM" || model == "SPGM" || model == "TPGM"){
		## implement the speed.pois="fast" in mp.generator
		lambda = 2
		lambda.c = 0.5
		
		A = adj2A(B,type=graph.type)
		sigma=lambda*B
		ltri.sigma=sigma[lower.tri(sigma)]
		nonzero.sigma = ltri.sigma[which(ltri.sigma !=0 )]
		Y.lambda = c(rep(lambda,nrow(sigma)), nonzero.sigma)
		
		Y = rmpois(n,Y.lambda)
		X =A%*%Y
		# add the labmda.c to all the nodes. 
		X = X + rmpois(ncol(X), rep(lambda.c,nrow(X)))
		mydata$X=X
		##mydata$A=A
		##mydata$sigma=sigma
	}
	#
	if(model == "GGM"){
		# Simmuate a Gaussian multivariate data matrix
		u = 0.1
		v = 0.3
		diag(B) = 0
		omega = B * v
		diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
		sigma = cov2cor(solve(omega))
		omega = solve(sigma)
		X = MASS::mvrnorm(n, rep(0, p), sigma)
		
		mydata$X <- t(X)
		##mydata$A=NULL
		##mydata$sigma=sigma
	}
	#
	if(model == "ISM"){
		theta = 0.4
		maxit = 1000
		
		X <- matrix(rbinom(n*p, 1, 0.6), nrow=n, ncol=p)
		X[X==0] <- -1
		# Changing sign of the parts of adjacency matrix 
		####dims = 1:size(x,2);Ising0.
		Theta <- -theta*B
		Theta[,1:(p/2)] <- -Theta[,1:(p/2)]
		Iterat <- array(0, dim=c(n,p,maxit))
		t = 1
		# RUN GIBBS SAMPLER
		while (t < maxit){
			t = t+1
			# Loop over dimensions  
			for (iD in 1:p){
				odds <- exp(2*X[, -iD]%*%Theta[-iD, iD])
				prob <- odds /(1+odds)
				X[, iD] <- rbinom(length(prob), 1, prob)*2 -1  
			}
			Iterat[,,t] <- X
		}
		xDat = (X+1)/2
		mydata$X <- t(xDat)
		##mydata$sigma=NULL
	}
	#
	if(model == "PGM"){
		mu = 0.5
		sigma = sqrt(0.01)
		
    mu = 0.1
    signam = 0.01
    
		Theta = matrix(rnorm(p*p,mu,sigma),nrow = p)
		Theta = Theta * B
		X = PGMSim(n,p, alpha = rep(0,p), Theta=Theta, maxit=1000)
		X = PGMSim(n,p, alpha = rep(0,p), Theta=Theta, maxit=130)
		
		mydata$X=X
		##mydata$A=NULL
		##mydata$sigma=NULL
	}

	return(mydata)
}
