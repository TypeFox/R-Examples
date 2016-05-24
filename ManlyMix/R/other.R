Gmat <- function(p){

    M <- p * (p + 1) / 2

    G <- matrix(c(rep(NA, M * p^2)), ncol = M)

    n <- 1

    Z <- rep(0, M)

    for (a in 1:p){

        for (b in 1:p){

            Zn <- Z

            i1 <- max(a, b)            
            i2 <- min(a, b)

            ind <- M - (p - i2 + 1) * (p - i2 + 2) / 2 + i1 - i2 + 1

            Zn[ind] <- 1

            G[n,] <- Zn

            n <- n + 1
            
        }

    }

    return(G)

}



Manly.contour <- function(X, var1 = 1, var2 = 2, model = NULL, x.slice = 100, y.slice = 100, x.mar = 1, y.mar = 1, col = "lightgrey", ...){
	n <- dim(X)[1]
	p <- dim(X)[2]
 	if(p != 2) warning("Only two variables are used for the contour plot...\n")

	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")


	x.seq <- seq(min(X[,var1]) - x.mar, max(X[,var1]) + x.mar, length = x.slice)
	y.seq <- seq(min(X[,var2]) - y.mar, max(X[,var2]) + y.mar, length = y.slice)
	z <- matrix(0, length(x.seq), length(y.seq))

	for (i in 1:length(x.seq)){
		for (j in 1:length(y.seq)){
			z[i,j] <- Manly.mix(c(x.seq[i], y.seq[j]), la = model$la[,c(var1, var2)], tau = model$tau, Mu = model$Mu[,c(var1, var2)], S = model$S[c(var1, var2),c(var1, var2),])
		}
	}	

	contour(x.seq, y.seq, z, col = col, ...)
	points(X[,c(var1, var2)], col=model$id, cex = 0.5, pch = 19)

}

Manly.mix <- function(X, la = NULL, tau = NULL, Mu = NULL, S = NULL){


	if(is.null(dim(X))){
		n <- 1
		p <- dim(la)[2]

	}
	else{
		n <- dim(X)[1]
		p <- dim(X)[2]
		
	}
	K <- length(tau)

	if( is.null(Mu) || is.null(la) || is.null(tau) || is.null(S)) stop("Must provide the parameters to calculate mixture density from...\n")

	equal.K <- c(dim(la)[1], dim(Mu)[1], dim(S)[3])
	equal.p <- c(dim(la)[2], dim(Mu)[2], dim(S)[1], dim(S)[2])


	if(K < 1) stop("Wrong number of mixture components K...\n")
	if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3])) stop("Inconsistent number of mixture components K...\n")
	if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionaltiy p...\n")


	x1 <- as.vector(t(X))
	la1 <- as.vector(t(la))
	Mu1 <- as.vector(t(Mu))
	S1 <- as.vector(S)


	misc_int <- c(p, n, K)
	mix1 <- rep(0, n)


	result <- .C("run_Manly_mix", x1 = as.double(x1), misc_int = as.integer(misc_int), tau = as.double(tau), Mu1 = as.double(Mu1), S1 = as.double(S1), la1 = as.double(la1), mix = as.double(mix1), PACKAGE = "ManlyMix")	
	
	return(result$mix)



}



Manly.var <- function(X, model = NULL){

	n <- dim(X)[1]
	p <- dim(X)[2]
	K <- length(model$tau)
	
	Inf.mat <- NULL

	tau <- model$tau
	gamma <- model$gamma
	Mu <- model$Mu
	S <- model$S
	la <- model$la

	x1 <- as.vector(t(X))


	misc_int <- c(p, n, K)
	
	Y <- array(NA, c(n, p, K))

	for(k in 1:K){

		y1 <- rep(0, n*p)
		result <- .C("run_Manly_transX", x1 = as.double(x1), misc_int = as.integer(misc_int), la1 = as.double(la[k,]), y1 = as.double(y1), PACKAGE = "ManlyMix")
		Y[,,k] <- matrix(result$y1, nrow = n, byrow = TRUE)

	}


	for (i in 1: n){

		grad1 <- NULL

		for (k in 1:(K-1)){

			grad1 <- c(grad1, gamma[i,k] / tau[k] - gamma[i,K] / tau[K])

		}


		grad2 <- NULL


		for (k in 1:K){
			grad2 <- c(grad2, gamma[i,k] * solve(S[,,k]) %*% (Y[i,,k] - Mu[k,]))

		}


		grad3 <- NULL

		I <- matrix(0,p,p)
		diag(I) <- rep(1,p)

		for (k in 1:K){

			grad3 <- c(grad3, t(Gmat(p)) %*% as.vector(gamma[i,k] / 2 * solve(S[,,k]) %*% ((Y[i,,k] - Mu[k,]) %*% t(Y[i,,k] - Mu[k,]) %*% solve(S[,,k]) - I)))


		}

		grad4 <- NULL


		for (k in 1:K){

			Dk <- matrix(0,p,p)
			for (j in 1:p){
				if(la[k,j] != 0){
					Dk[j,j] <- (1 + exp(la[k,j] * X[i,j]) * (X[i,j] * la[k,j] - 1)) / la[k,j]^2
				}
			}

			vect.temp <- -gamma[i,k] * Dk %*% solve(S[,,k]) %*% (Y[i,,k] - Mu[k,]) + gamma[i,k] * X[i,]
			index <- la[k,] == 0
			grad4 <- c(grad4, vect.temp[!index])


		}
		q <-  c(grad1, grad2, grad3, grad4)
	
		
		Inf.mat <- cbind(Inf.mat, q)

	}
	Inf.mat <- Inf.mat %*% t(Inf.mat)

	V <- solve(Inf.mat)

	return(V)
}




Manly.sim <- function(n, la, tau, Mu, S){

	if(n < 1)  stop("Wrong sample size n...\n")
	if(sum(tau) != 1)  stop("Mixing proportions should sum up to 1...\n")

    	K <- length(tau)
    	p <- dim(Mu)[2]

	equal.K <- c(dim(la)[1], dim(Mu)[1], dim(S)[3])
	equal.p <- c(dim(la)[2], dim(S)[1], dim(S)[2])

	if (K < 1) stop("Wrong number of mixture components K...\n")
	if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3])) stop("Inconsistent number of mixture components K...\n")
	if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3])) stop("Inconsistent number of dimensionality p...\n")

	Y <- NULL
	Y <- matrix(rnorm(n * p), n, p)
        Nk <- rep(1, K) + drop(rmultinom(1, n - K, tau))
   	
    	id <- NULL
    	for (k in 1:K) {id <- c(id, rep(k, Nk[k]))}
	for (k in 1:K) {
	   	eS <- eigen(S[,,k], symmetric = TRUE)
   		ev <- eS$values	
		temp <- eS$vectors %*% diag(sqrt(ev), p) %*% t(Y[id == k, ]) + Mu[k,]
        	Y[id == k, ] <- t(temp)

	}


	X <- NULL	
	for (k in 1:K){
		ind <- id == k
		Z <- sweep(Y[ind,], 2, STATS = la[k,], FUN = "*")
		Z <- log(Z + 1)
		Z <- sweep(Z, 2, STATS = la[k,], FUN = "/")
		X <- rbind(X, Z)
	}

	return(list(X = X, id = id))
}	




