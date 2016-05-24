## This is the R/Splus programs for 
## homogenuous dispersion simplex marginal model (song and Tan, 2000) 
## and heterogeneuous dispersion simplex marginal model (song, Qiu and Tan, 2004)

## note: sigma represent the dispersion parameter sigma^2

#### Generic functions:
DLOW <- function(x) {
	r <- dim(x)[1]
	for(j in 1:r) {
		for(k in 1:r) {
			if(j < k) {
				x[k, j] <- x[k, j]
			}
			else x[k, j] <-  - Inf
		}
	}
	bx <- as.vector(x)
	bx <- bx[bx >  - Inf]
	return(bx)
}

norm <- function(x) {
# x is a vector, to get the type 2-norm of x
     sqrt(sum(x^2))    
}

halfEd2 <- function(mu, sigma) {
# variance of score function
    (3 * sigma)/(mu * (1 - mu)) + 1/Vfun(mu)
}

#### GEE for heterogeneity dispersion simplex marginal model (Song, Qiu and Tan, 2004)

simgee2 <- function(y, X, Z, T, id, link = 1, type = 0, beta, Sigma, alpha, maxiter, tol) {
## This program is used to produce the estimates of parameters for  
## AR(1) correlation structure by GEE method, where alpha is estimated
## by an estimating equation, for simplex marginal model with heterogeneity-dispersion
## model: y_ij~Simplex(mu_ij,sigma_ij), 
##       logit(mu_ij)=beta0+beta1*x_1ij+beta2*x_2ij
##       sigma_ij = sigma0 + sigma1*z_1ij + sigma2*z_2ij
## Attached functions: gf, gfd, hf, tf, wwf, dd, INV, Vfun, uu, halfEd2, DLOW
## sort the vector
   	Order <- order(id)
   	y <- y[Order]
   	X <- X[Order, ]
   	Z <- Z[Order, ]
   	T <- T[Order]
	m <- length(levels(id))
## make an n_i vector:
	ID <- rep(0, m)
	for(i in 1:m) {
		ID[i] <- sum(id == (levels(id)[i]))
	}
	ID <- cumsum(ID)
   	N <- length(y)
	q <- dim(Z)[2]
	p <- dim(X)[2]	
## specify the link function
    if (link==1){
		gf <- function(x) gf1(x)
		hf <- function(x) hf1(x)
		gfd <- function(x) gfd1(x)
	}
	else if (link==2) {
		gf <- function(x) gf2(x)
		hf <- function(x) hf2(x)
		gfd <- function(x) gfd2(x)
	}
	else if (link==3) {
		gf <- function(x) gf3(x)
		hf <- function(x) hf3(x)
   	   	gfd <- function(x) gfd3(x)
    }
   	else if (link==4) {
   	   	gf <- function(x) gf4(x)
   	   	hf <- function(x) hf4(x)
   	   	gfd <- function(x) gfd4(x)
   	}
   	else {cat("Link function is not defined")}
## initial values of beta and mu:
   	sim.glm.out <- simglm1(y, X, link = link, beta = NULL, maxit = 200, tol = 0.1^6)
   	if (length(beta) != p)
		beta <- as.vector(sim.glm.out$fixef[,1])
	mu <- as.vector(hf(X %*% beta))
## initial values of sigma:
	ddd <- dd(y, mu)
	options(warn = -1)
   	out <- try(glm(ddd ~ Z[,-1], family = Gamma(link = log)), TRUE)
   	if(inherits(out, "try-error"))
   	   	out <- glm(ddd ~ Z[,-1], family = Gamma(link = log), start = c(mean(log(ddd)), 
			rep(0, q-1)), maxit = 5000)
	if (length(Sigma) != q)		
		Sigma <- as.numeric(out$coef)
	sigma <- as.vector(exp(Z %*% Sigma))
## initial values of alpha
	rr <- uu(y, mu) / sqrt(sigma * halfEd2(mu, sigma))
	if (is.null(alpha)){
		alpha <- -0.5
	}
## Newton-Scoring algorithm for theta:
	theta <- c(beta, Sigma, alpha)
#	cat("initial theta-Sigma is", round(c(theta, Sigma), 4), "\n")
	iter <- 0
	dif <- 1
	theta.old <- theta
	while(dif > tol && iter < maxiter) {
		iter <- iter + 1
   	   	result <- .C("Update", as.double(beta), as.double(Sigma), as.double(alpha), as.double(mu), 
   	   	   	as.double(sigma), as.double(y), as.double(t(X)), as.double(t(Z)), as.double(T), as.integer(ID), 
			as.integer(N), as.integer(m), as.integer(p), as.integer(q), as.integer(link), as.integer(type))
   	   	beta <- result[[1]]
   	   	Sigma <- result[[2]]
   	   	alpha <- result[[3]]
		if (alpha== -1000){
			return(NULL)
		}
   	   	mu <- result[[4]]
   	   	sigma <- result[[5]]
   	   	theta <- c(beta, Sigma, alpha)
   	   	diff <- as.vector(theta - theta.old)
   	   	dif <- norm(diff) / norm(theta)
   	   	theta.old <- theta
   	}
	U <- uu(y, mu)
	U1 <- sqrt(sigma * halfEd2(mu, sigma))
	DD <- dd(y, mu) - sigma
	V11 <- matrix(0, p, p)
	V12 <- matrix(0, p, q)
	V13 <- matrix(0, p, 1)
	V22 <- matrix(0, q, q)
	V23 <- matrix(0, q, 1)
	V33 <- matrix(0, 1, 1)
	S1 <- matrix(0, p, p)
	S3 <- 0
	for(i in 1:m) {
		if (i == 1)
			a <- 1
		else
			a <- ID[i-1] + 1
		b <- ID[i]
		n <- b - a + 1
		w <- 1.0 / gfd(mu[a:b])
		D <- t(t(X[a:b,  ]) %*% diag(w))
		ww <- as.vector(3 * sigma[a:b] * pp(mu[a:b])^2 + 1)
		A <- diag(ww)
		www <- 1/sqrt(sigma[a:b] * Vfun(mu[a:b]) * ww)
		DVf <- diag(Vfun(mu[a:b]))
		dR <- matrix(0, n, n)
		Eta <- diag(n)
   	   	if (type == 0){
   	   	   	for(j in 1:(n - 1)) {
				for(k in (j + 1):n) {
					Eta[k, j] <- exp(alpha * abs(T[a + k - 1] - 
						T[a + j - 1]))
					dR[k, j] <- abs(T[a + k - 1] - T[a + j - 
						1]) * Eta[k, j]
				}
			}
		}
		else{
			Eta <- matrix(exp(alpha) / 2, n, n)
			diag(Eta) <- 1
			dR <- matrix(exp(alpha) / 2, n, n)
			diag(dR) <- 0
		}
		dR <- dR + t(dR)
		cc <- DLOW(dR)
		Eta <- Eta + t(Eta)
		Ures <- outer(U[a:b] / U1[a:b], U[a:b] / U1[a:b])
		rr <- DLOW(Ures) - DLOW(Eta)
		V <- INV(Eta + t(Eta) - diag(n))
		VV <- diag(www) %*% V %*% diag(www)
		S1 <- S1 + t(D) %*% A %*% VV %*% A %*% D
		S3 <- S3 - as.vector(DLOW(dR) %*% DLOW(dR))
		UU <- outer(U[a:b], U[a:b])
		Ud <- outer(U[a:b], DD[a:b])
		Ur <- outer(U[a:b], rr)
		Czz <- DVf %*% UU %*% DVf
		Czd <- DVf %*% Ud
		Czr <- DVf %*% Ur
		Cdd <- outer(DD[a:b], DD[a:b])
		Cdr <- outer(DD[a:b], rr)
		Crr <- outer(rr, rr)
		ZZ <- diag(1.0 / sigma[a:b]) %*% Z[a:b, ] / 2
		V11 <- V11 + t(D) %*% A %*% VV %*% Czz %*% VV %*% A %*% D
		V12 <- V12 + t(D) %*% A %*% VV %*% Czd %*% ZZ
		V13 <- V13 + t(D) %*% A %*% VV %*% Czr %*% cc
		V22 <- V22 + t(ZZ) %*% Cdd %*% ZZ
		V23 <- V23 + t(ZZ) %*% Cdr %*% cc
		V33 <- V33 + t(cc) %*% Crr %*% cc
	}
	S2 <- t(Z) %*% Z / 2
	S11 <- INV(S1)
	S22 <- INV(S2)
	S33 <- 1/S3
	SS <- matrix(0, (p + q + 1), (p + q + 1))
	SS[1:p, 1:p] <- S11
	SS[(p + 1):(p + q), (p + 1):(p + q)] <- S22
	SS[(p + q + 1), (p + q + 1)] <- S33
	V21 <- t(V12)
	V31 <- t(V13)
	V32 <- t(V23)
	VVV <- matrix(0, (p + q + 1), (p + q + 1))
	VVV[1:p, 1:p] <- V11
	VVV[1:p, (p + 1):(p + q)] <- V12
	VVV[1:p, (p + q + 1)] <- V13
	VVV[(p + 1):(p + q), 1:p] <- V21
	VVV[(p + 1):(p + q), (p + 1):(p + q)] <- V22
	VVV[(p + 1):(p + q), (p + q + 1)] <- V23
	VVV[(p + q + 1), 1:p] <- V31
	VVV[(p + q + 1), (p + 1):(p + q)] <- V32
	VVV[(p + q + 1), (p + q + 1)] <- V33
	invGinf <- SS %*% VVV %*% SS
	std <- sqrt(diag(invGinf))
	stdfixef <- std[1:p]
   	stdscor <- U / U1
	stddispar <- std[(p + 1):(p + q)]
	stdalpha <- std[(p + q + 1)]
	stdcorpar <- std[(p + q + 1)] * exp(alpha)	
	eee <- as.vector((y - mu)/sqrt(var.sim(mu, sigma)))
	ee <- as.vector((y - mu)/sqrt(mu * (1 - mu)))
	ss <- sigma * ((3 * sigma)/(mu * (1 - mu)) + 1/Vfun(mu))
	sss <- as.vector(gf(mu) + uu(y, mu)/sqrt(ss))
	devi <- dd(y,mu)/sigma
	pred <- as.vector(X %*% beta)
	options(warn = 0)
	if (iter == maxiter)
		warning("step size truncated due to divergence")
	return(list(fixef = cbind(beta, stdfixef), dispar = cbind(Sigma, stddispar), 
		Dispersion = sigma, autocor = c(alpha, stdalpha, exp(alpha), stdcorpar),
		appstdPerr = ee, stdPerr = eee, meanmu = mu, adjvar = sss, stdscor = stdscor, 
   	   	deviance = devi, predict = pred, iter = iter))
}

#### GEE for homogeneous dispersion simplex marginal model (Song and Tan, 2000)

simgee1 <- function(y, X, T, id, link = 1, type = 0, beta, alpha, maxiter, tol) {
##This program is used to produce the estimates of parameters for the case 
##of AR(1) correlation structure for GEE method, where  alpha is estimated
##by an estimating equation. Assume homogeneity in dispersion parameter.
##model: y_ij~Simplex(mu_ij,sigma), 
##       logit(mu_ij)=beta0+beta1*x_1ij+beta2*x_2ij
##Attached functions: gf, gfd, hf, tf, wwf, dd, INV, Vfun, uu, halfEd2, DLOW
## sort the vector
   	Order <- order(id)
   	y <- y[Order]
   	X <- X[Order, ]
   	T <- T[Order]
	m <- length(levels(id))
## make an n_i vector:
	ID <- rep(0, m)
	for(i in 1:m) {
		ID[i] <- sum(id == levels(id)[i])
	}
	ID <- cumsum(ID)
	N <- length(y)
	Y <- cbind(rep(1, N))
	p <- dim(as.matrix(X))[2]  
## specify the link function
    if (link==1){
		gf <- function(x) gf1(x)
		hf <- function(x) hf1(x)
		gfd <- function(x) gfd1(x)
		}
	else if (link==2) {
		gf <- function(x) gf2(x)
		hf <- function(x) hf2(x)
		gfd <- function(x) gfd2(x)
		}
	else if (link==3) {
		gf <- function(x) gf3(x)
		hf <- function(x) hf3(x)
   	   	gfd <- function(x) gfd3(x)
   	   	}
   	else if (link==4) {
   	   	gf <- function(x) gf4(x)
   	   	hf <- function(x) hf4(x)
   	   	gfd <- function(x) gfd4(x)
   	   	}
   	else {cat("Link function is not defined")}
## initial values of beta and mu:
   	sim.glm.out <- simglm1(y, X, link = link, beta = NULL, maxit = 200, tol = 0.1^6)
	if (length(beta) != p)
		beta <- as.vector(sim.glm.out$fixef[,1])
	mu <- as.vector(hf(X %*% beta))
	sigma <- rep(sum(dd(y, mu)) / (N - p), N)
## initial values of alpha:
	rr <- uu(y, mu) / sqrt(sigma * halfEd2(mu, sigma))
	if (is.null(alpha)){
#		alpha <- preCor(T, rr, ID)
#		if (alpha >= 0||is.na(alpha))
		alpha <- -0.5
	}
## Newton-Scoring algorithm for theta:
	theta <- c(beta, alpha)
#	cat("initial theta-Sigma is", round(c(theta, Sigma), 4), "\n")
	iter <- 0
	dif <- 1
	theta.old <- theta
	while(dif > tol && iter < maxiter) {
		iter <- iter + 1
   	   	result <- .C("Update_homo", as.double(beta), as.double(alpha), as.double(mu), 
			as.double(sigma), as.double(y), as.double(t(X)), as.double(T), as.integer(ID), 
			as.integer(N), as.integer(m), as.integer(p), as.integer(link), as.integer(type))
   	   	beta <- result[[1]]
   	   	alpha <- result[[2]]
		if (alpha== -1000){
			return(NULL)
		}
   	   	mu <- result[[3]]
   	   	sigma <- result[[4]]
   	   	theta <- c(beta, alpha)
   	   	diff <- as.vector(theta - theta.old)
   	   	dif <- norm(diff) / norm(theta)
   	   	theta.old <- theta
	}
	U <- uu(y, mu)
	U1 <- sqrt(sigma * halfEd2(mu, sigma))
	V11 <- matrix(0, p, p)
	V13 <- matrix(0, p, 1)
	V33 <- matrix(0, 1, 1)
	S1 <- matrix(0, p, p)
	S3 <- 0
	for(i in 1:m) {
		if (i == 1)
			a <- 1
		else
			a <- ID[i-1] + 1
		b <- ID[i]
		n <- b - a + 1
		w <- 1.0 / gfd(mu[a:b])
		D <- t(t(X[a:b,  ]) %*% diag(w))
		ww <- as.vector(3 * sigma[a:b] * pp(mu[a:b])^2 + 1)
		A <- diag(ww)
		www <- 1/sqrt(sigma[a:b] * Vfun(mu[a:b]) * ww)
		DVf <- diag(Vfun(mu[a:b]))
		dR <- matrix(0, n, n)
		Eta <- diag(n)
		if (type == 0){
			for(j in 1:(n - 1)) {
				for(k in (j + 1):n) {
					Eta[k, j] <- exp(alpha * abs(T[a + k - 1] - 
						T[a + j - 1]))
					dR[k, j] <- abs(T[a + k - 1] - T[a + j - 
						1]) * Eta[k, j]
				}
			}
		}
		else{
			Eta <- matrix(exp(alpha) / 2, n, n)
			diag(Eta) <- 1
			dR <- matrix(exp(alpha) / 2, n, n)
			diag(dR) <- 0
		}
		dR <- dR + t(dR)
		cc <- DLOW(dR)
		Eta <- Eta + t(Eta)
		Ures <- outer(U[a:b] / U1[a:b], U[a:b] / U1[a:b])
		rr <- DLOW(Ures) - DLOW(Eta)
		V <- INV(Eta + t(Eta) - diag(n))
		VV <- diag(www) %*% V %*% diag(www)
		S1 <- S1 + t(D) %*% A %*% VV %*% A %*% D
		S3 <- S3 - as.vector(DLOW(dR) %*% DLOW(dR))
		UU <- outer(U[a:b], U[a:b])
		Ur <- outer(U[a:b], rr)
		Czz <- DVf %*% UU %*% DVf
		Czr <- DVf %*% Ur
		Crr <- outer(rr, rr)
		V11 <- V11 + t(D) %*% A %*% VV %*% Czz %*% VV %*% A %*% D
		V13 <- V13 + t(D) %*% A %*% VV %*% Czr %*% cc
		V33 <- V33 + t(cc) %*% Crr %*% cc
	}
	S11 <- INV(S1)
	S33 <- 1/S3
	SS <- matrix(0, (p + 1), (p + 1))
	SS[1:p, 1:p] <- S11
	SS[(p + 1), (p + 1)] <- S33
	V31 <- t(V13)
	VVV <- matrix(0, (p + 1), (p + 1))
	VVV[1:p, 1:p] <- V11
	VVV[1:p, (p + 1)] <- V13
	VVV[(p + 1), 1:p] <- V31
	VVV[(p + 1), (p + 1)] <- V33
	invGinf <- SS %*% VVV %*% SS
	std <- sqrt(diag(invGinf))
	stdfixef <- std[1:p]
   	stdscor <- U / U1
	stdalpha <- std[(p + 1)]
	stdcorpar <- std[(p + 1)] * exp(alpha)
    eee <- as.vector((y - mu)/sqrt(var.sim(mu, sigma)))
	ee <- as.vector((y - mu)/sqrt(mu * (1 - mu)))
	ss <- sigma * ((3 * sigma)/(mu * (1 - mu)) + 1/Vfun(mu))
	sss <- as.vector(gf(mu) + uu(y, mu)/sqrt(ss))
	devi <- dd(y, mu) / sigma[1]
	pred <- as.vector(X %*% beta)
	if (iter == maxiter)
		warning("step size truncated due to divergence")
	return(list(fixef = cbind(beta, stdfixef), Dispersion = sigma[1], autocor = c(alpha, 
		stdalpha, exp(alpha), stdcorpar), appstdPerr = ee, stdPerr = eee, stdscor = stdscor, 
   	   	meanmu = mu, adjvar = sss, deviance = devi, predict = pred, iter = iter))
}

cormm <- function(id, T, stdscor) {
	m <-length(levels(id))
	mm <- matrix(0, m, max(T))
	for(i in 1:m) {
		a <- T[id == levels(id)[i]]
		b <- stdscor[id == levels(id)[i]]
		mm[i, a] <- b
	}
	return(mm)
}

plotmm <- function(mm, r, ...){
   	m <- dim(mm)[1]
   	T <- dim(mm)[2]
   	x <- integer(0)
   	y <- integer(0)
   	for (i in 1:m){
   	   	s <- ((mm[i,1:(T-r)] != 0) & (mm[i,(r+1):T] != 0))
   	   	x <- c(x, mm[i,1:(T-r)][s])
   	   	y <- c(y, mm[i,(r+1):T][s])
   	}
	bound <- max(abs(x), abs(y))
	plot(x, y, ...)
}

preCor <- function(T, rr, ID){
	m <- length(ID)
	su_cor <- 0
	s <- 0
	for (i in 1:m){
		if (i == 1)
			a <- 1
		else
			a <- ID[i-1]+1
		b <- ID[i]
		corre <- outer(rr[a:b], rr[a:b])
		for (k in a:(b-1)){
			for (j in (k+1):b){
				if (abs(T[k] - T[j]) == 1){
					su_cor <- su_cor + corre[k-a+1, j-a+1] 
					s <- s + 1
				}
			}
		}
	}
	if (s == 0)
		return (-0.5)
	else
		return (log(su_cor / s))
}