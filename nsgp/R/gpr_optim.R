

# regular 3-parameter optimisation function
.optim_stat <- function(params, x,y, x.targets, noise.targets, alpha=1, beta=0) {
	alpha + beta # annoying R, do some fake stuff with alpha and beta
	
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	n <- length(x)
	y <- y - mean(y)
	
	# compute kernel
	K.y <- .RBF(x, x, sigma.f, l) + diag(noise.targets[match(x, x.targets)])^2 # add noise term
	
	return (as.numeric(-0.5* t(y) %*% solve(K.y) %*% y -0.5* log(det(K.y)) - 0.5*n*log(2*pi)))
}


# dynamic 5-parameter optimisation function
.optim_dyn <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0) {
	alpha + beta # annoying R, do some fake stuff with alpha and beta
	
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	lmin    <- params[4]
	c       <- params[5]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	
	n <- length(x)
	y <- y - mean(y)
	
	# compute kernel
	K.y <- .nsRBF(x, x, sigma.f, l, lmin, c) + diag(noise.targets[match(x, x.targets)])^2 # add noise term
	
	return (as.numeric(-0.5* t(y) %*% solve(K.y) %*% y -0.5* log(det(K.y)) - 0.5*n*log(2*pi)))
}

# dynamic 5-parameter optimisation function for greedy data fit
.optim_fit_dyn <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0, lambda=0.5) {
	alpha + beta # annoying R, do some fake stuff with alpha and beta
	
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	lmin    <- params[4]
	c       <- params[5]
	
	if (all(is.na(noise.targets))) {
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	
	n <- length(x)
	y <- y - mean(y)
	
	# compute kernel
	K.y <- .nsRBF(x, x, sigma.f, l, lmin, c) + diag(noise.targets[match(x, x.targets)])^2 # add noise term
	
	# double the effect of data fit, effectively \lambda=0.5
	return (as.numeric(-0.5*t(y)%*%solve(K.y)%*%y -0.5*lambda*log(det(K.y)) -0.5*n*log(2*pi)))
}


.optim_exp <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0) {
	alpha + beta # annoying R, do some fake stuff with alpha and beta
		
	if (length(params) == 4) {
		params = c(params[1], 1, params[2:4])
	}
	
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	lmin    <- params[4]
	c       <- params[5]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	
	n2 <- length(x.targets)
	y <- y - mean(y)
	
	# compute kernels
	K = .nsRBF.blocks(x, x.targets, sigma.f, l, lmin, c)
	K.tt = K$tt
	K.ts = K$ts
	K.st = K$st
	K.ss = K$ss
	K.y <- K.tt + diag(noise.targets[ match(x, x.targets) ])^2 # add noise variance
	K.y.inv <- solve(K.y)
	
	# Calculate the models
	f.mean <- K.st%*%K.y.inv%*%y
	f.cov <- K.ss - K.st%*%K.y.inv%*%K.ts
	K.m <- f.cov + K.ss + 2*diag(noise.targets)^2

	
	val <- as.numeric(-0.5*t(f.mean)%*%solve(K.m)%*%f.mean -0.5*log(det(K.m)) -0.5*n2*log(2*pi))

#	cat('optim_exp:  ')
	loss = -0.5*t(f.mean)%*%solve(K.m)%*%f.mean
	reg = -0.5*log(det(K.m))
	const = -0.5*n2*log(2*pi)
#	cat("EMLL=", val, ' loss=', loss , ' reg=', reg, ' const=', const)
#	cat(' params=', params, '\n')
	
	return (val)
}

.optim_exp_sample <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0) {
	alpha + beta # annoying R, do some fake stuff with alpha and beta
	
	if (length(params) == 4) {
		params = c(params[1], 1, params[2:4])
	}
	
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	lmin    <- params[4]
	c       <- params[5]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	
	n2 <- length(x.targets)
	y <- y - mean(y)
	
	# compute kernels
	K = .nsRBF.blocks(x, x.targets, sigma.f, l, lmin, c)
	K.tt = K$tt
	K.ts = K$ts
	K.st = K$st
	K.ss = K$ss
	K.y <- K.tt + diag(noise.targets[ match(x, x.targets) ])^2 # add stdev noise term
	K.y.inv <- solve(K.y)
	
	# Calculate the models
	f.mean <- K.st%*%K.y.inv%*%y
	#  f.cov <- K.ss - K.st%*%K.y.inv%*%K.ts
	K.m <- K.ss + diag(noise.targets)^2
	
	val <- as.numeric(-0.5*t(f.mean)%*%solve(K.m)%*%f.mean -0.5*log(det(K.m)) -0.5*n2*log(2*pi))
	return (val)
}


# gradients of the 3 hyperparameters in the standard case
.optim_stat_gradients <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0) {
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	
	n <- length(y)
	y <- y - mean(y)
	
	# compute kernels
	K.tt <- .RBF(x, x, sigma.f, l)
	K.y <- K.tt + diag(noise.targets[match(x,x.targets)])^2 # add stdev noise term
	K.y.inv <- solve(K.y)
	
	gKf <- 2*K.tt/sigma.f
	gKl <- mat.or.vec(n,n)
	for (i in 1:n) {
		for (j in 1:n) {
			t1 <- x[i]
			t2 <- x[j]
			gKl[i,j] <- 2*(t1-t2)^2/l^3 * K.tt[i,j]
		}
	}
	
#	gKn = diag(2*noise.targets[match(x,x.targets)]^2)
	
	alphas <- K.y.inv %*% y
	
	gf    <- 0.5* sum(diag((alphas %*% t(alphas) - K.y.inv) %*% gKf))
	gn    <- 0.5* sum(diag((alphas %*% t(alphas) - K.y.inv) %*% (2*sigma.n*diag(n)))) + (alpha-1)/sigma.n - beta
	gl    <- 0.5* sum(diag((alphas %*% t(alphas) - K.y.inv) %*% gKl))
	
	return (c(gf,gn,gl))
}

# gradients of the 5 hyperparameters in the dynamic case
.optim_dyn_gradients <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0) {
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	lmin    <- params[4]
	c       <- params[5]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	
	n <- length(y)
	y <- y - mean(y)
	
	# compute kernels
	O.tt    = diag(noise.targets[match(x, x.targets)])^2
	K.tt    = .nsRBF(x, x, sigma.f, l, lmin, c) # data cov
	K.y     = K.tt + O.tt # add stdev noise term
	K.y.inv = solve(K.y)
	
	# shortcuts
	tt = gradientfuncs(x, x, l, lmin, c)
	
	# gradient matrices for kernel K.y
	g.Ky.f    = 2 / sigma.f * K.tt
	g.Ky.n    = 2 * sqrt(O.tt)
	g.Ky.l    = -2 * tt$A * tt$B * K.tt
	g.Ky.lmin = -2 * tt$A * tt$C * K.tt
	g.Ky.c    = -2 * tt$A * tt$D * K.tt
	
	alphas <- K.y.inv %*% y
	
	gf    <- 0.5*sum(diag((alphas %*% t(alphas) - K.y.inv) %*% g.Ky.f))
	gn    <- 0.5*sum(diag((alphas %*% t(alphas) - K.y.inv) %*% g.Ky.n)) + (alpha-1)/sigma.n - beta
	gl    <- 0.5*sum(diag((alphas %*% t(alphas) - K.y.inv) %*% g.Ky.l))
	glmin <- 0.5*sum(diag((alphas %*% t(alphas) - K.y.inv) %*% g.Ky.lmin))
	gc    <- 0.5*sum(diag((alphas %*% t(alphas) - K.y.inv) %*% g.Ky.c))
	
	return (c(gf,gn,gl,glmin,gc))
}

# gradients of the 5 hyperparameters in the dynamic case
# the 'lambda' parameter gives additional weight to either fit or regularization
.optim_fit_dyn_gradients <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0, lambda=0.5) {
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	lmin    <- params[4]
	c       <- params[5]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	} 
	
	n <- length(y)
	y <- y - mean(y)
	
	# compute kernels
	K.tt <- .nsRBF(x, x, sigma.f, l, lmin, c) # data cov
	K.y <- K.tt + diag(noise.targets[ match(x, x.targets) ])^2 # add stdev noise term
	K.y.inv <- solve(K.y)
	
	# shortcuts
	bexp2 <- function(x) { .bexp(x, l, lmin, c) }
	s1 <- function(t1,t2) {    return (t1/bexp2(t1) - t2/bexp2(t2))  }
	s2 <- function(t1,t2) {    return ((t2-t2*exp(-c*t2))/bexp2(t2)^2 - (t1-t1*exp(-c*t1))/bexp2(t1)^2)  }
	s3 <- function(t1,t2) {    return ((t2*exp(-c*t2))/bexp2(t2)^2 - (t1*exp(-c*t1))/bexp2(t1)^2)  }
	s4 <- function(t1,t2) {    return ( (l-lmin) * ( (t2^2*exp(-c*t2))/bexp2(t2)^2 - (t1^2*exp(-c*t1))/bexp2(t1)^2 ) )  }
	
	gKf <- 2/sigma.f * K.tt
	gKn <- 2*sigma.n*diag(n)
	gKl <- mat.or.vec(n,n)     # empty
	gKlmin <- mat.or.vec(n,n)  # empty
	gKc <- mat.or.vec(n,n)     # empty
	
	for (i in 1:n) {
		for (j in 1:n) {
			t1 <- x[i]
			t2 <- x[j]
			gKl[i,j]    <- -2*s2(t1,t2)*s1(t1,t2)*K.tt[i,j]
			gKlmin[i,j] <- -2*s3(t1,t2)*s1(t1,t2)*K.tt[i,j]
			gKc[i,j]    <- -2*s4(t1,t2)*s1(t1,t2)*K.tt[i,j]
		}
	}
	
	alphas <- K.y.inv %*% y
	
	gf    <- 0.5*t(alphas) %*% gKf %*% alphas    -0.5*lambda*sum(diag(K.y.inv %*% gKf))
	gn    <- 0.5*t(alphas) %*% gKn %*% alphas    -0.5*lambda*sum(diag(K.y.inv %*% gKn)) + (alpha-1)/sigma.n - beta
	gl    <- 0.5*t(alphas) %*% gKl %*% alphas    -0.5*lambda*sum(diag(K.y.inv %*% gKl))
	glmin <- 0.5*t(alphas) %*% gKlmin %*% alphas -0.5*lambda*sum(diag(K.y.inv %*% gKlmin))
	gc    <- 0.5*t(alphas) %*% gKc %*% alphas    -0.5*lambda*sum(diag(K.y.inv %*% gKc))
	
	return (c(gf,gn,gl,glmin,gc))
}

# helper function that returns the matrices A, B, C and D
gradientfuncs = function(x1, x2, l, lmin, c) {
	n = length(x1)
	n2 = length(x2)
	
	# curry a function with fixed parameters
	bexp2 <- function(x) { return(.bexp(x, l, lmin, c)) }
	
	s1 <- function(t1,t2) {    return (t1/bexp2(t1) - t2/bexp2(t2))  }
	s2 <- function(t1,t2) {    return ((t2-t2*exp(-c*t2))/bexp2(t2)^2 - (t1-t1*exp(-c*t1))/bexp2(t1)^2)  }
	s3 <- function(t1,t2) {    return ((t2*exp(-c*t2))/bexp2(t2)^2 - (t1*exp(-c*t1))/bexp2(t1)^2)  }
	s4 <- function(t1,t2) {    return ( (l-lmin) * ( (t2^2*exp(-c*t2))/bexp2(t2)^2 - (t1^2*exp(-c*t1))/bexp2(t1)^2) )  }
	
	A = mat.or.vec(n,n2)
	B = mat.or.vec(n,n2)
	C = mat.or.vec(n,n2)
	D = mat.or.vec(n,n2)
	
	for (i in 1:n) {
		for (j in 1:n2) {
			t1 = x1[i]
			t2 = x2[j]
			A[i,j] = s1(t1,t2)
			B[i,j] = s2(t1,t2)
			C[i,j] = s3(t1,t2)
			D[i,j] = s4(t1,t2)
		}
	}
	return(list(A=A,B=B,C=C,D=D))
}


.optim_exp_gradients <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0) {
#	par4 = FALSE
#	if (length(params) == 4) {
#		par4 = TRUE
#		params = c(params[1], 1, params[2:4])
#	}
	
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	lmin    <- params[4]
	c       <- params[5]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	
	n <- length(y)
	n2 <- length(x.targets)
	y <- y - mean(y)
	
	# compute kernels
	O.tt = diag(noise.targets[match(x, x.targets)])^2
	O.ss = diag(noise.targets)^2
	M = .nsRBF.blocks(x, x.targets, sigma.f, l, lmin, c)
	K.tt = M$tt
	K.ts = M$ts
	K.st = M$st
	K.ss = M$ss
	K.y     <- K.tt + O.tt # add noise variance term
	K.y.inv <- solve(K.y)
	f.mean  <- K.st%*%K.y.inv%*%y
	f.cov   <- K.ss - K.st%*%K.y.inv%*%K.ts
	K.m     <- f.cov + K.ss + 2*O.ss
	K.m.inv <- solve(K.m)
	
	# more stuff
	alphas <- K.m.inv %*% f.mean
	betas <- K.y.inv %*% K.ts
	
	# shortcuts
	E.ss = -2*K.st%*%betas + t(betas)%*%K.tt%*%betas + 2*K.ss
	
	tt = gradientfuncs(x, x, l, lmin, c)
	ts = gradientfuncs(x, x.targets, l, lmin, c)
	ss = gradientfuncs(x.targets, x.targets, l, lmin, c)
	
	# gradient matrices for kernel K.y
	g.Ky.f    = 2 / sigma.f * K.tt
	g.Ky.n    = 2 * sqrt(O.tt)
	g.Ky.l    = -2 * tt$A * tt$B * K.tt
	g.Ky.lmin = -2 * tt$A * tt$C * K.tt
	g.Ky.c    = -2 * tt$A * tt$D * K.tt
	
	# gradient matrices for kernel K.ts
	g.Kts.f    = 2 / sigma.f * K.ts
	g.Kts.n    = mat.or.vec(n,n2) # just zeros
	g.Kts.l    = -2 * ts$A * ts$B * K.ts
	g.Kts.lmin = -2 * ts$A * ts$C * K.ts
	g.Kts.c    = -2 * ts$A * ts$D * K.ts
	
	# gradient matrices for kernel K.m
	g.Km.f    <- 2 / sigma.f * E.ss
	g.Km.n    <- 2*t(betas)%*%sqrt(O.tt)%*%betas + 4*sqrt(O.ss)
	g.Km.l    <- -2 * ss$A * ss$B * E.ss
	g.Km.lmin <- -2 * ss$A * ss$C * E.ss
	g.Km.c    <- -2 * ss$A * ss$D * E.ss
	
	# gradients
	gf    <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% g.Km.f))    - t(alphas) %*% ((t(g.Kts.f)    - t(betas) %*% g.Ky.f) %*% (K.y.inv%*%y))
	gn    <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% g.Km.n))    - t(alphas) %*% ((t(g.Kts.n)    - t(betas) %*% g.Ky.n) %*% (K.y.inv%*%y)) + (alpha-1)/sigma.n - beta
	gl    <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% g.Km.l))    - t(alphas) %*% ((t(g.Kts.l)    - t(betas) %*% g.Ky.l) %*% (K.y.inv%*%y))
	glmin <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% g.Km.lmin)) - t(alphas) %*% ((t(g.Kts.lmin) - t(betas) %*% g.Ky.lmin) %*% (K.y.inv%*%y))
	gc    <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% g.Km.c))    - t(alphas) %*% ((t(g.Kts.c)    - t(betas) %*% g.Ky.c) %*% (K.y.inv%*%y))
	
#	cat('optim_exp_gradients:  ', c(gf,gn,gl,glmin,gc)/10000, '\n')
	
#	if (par4) { 
#		return (c(gf,gl,glmin,gc)) 
#	}
	
	return (c(gf,gn,gl,glmin,gc))
}

.optim_exp_sample_gradients <- function(params, x, y, x.targets, noise.targets, alpha=1, beta=0) {
	sigma.f <- params[1]
	sigma.n <- params[2]
	l       <- params[3]
	lmin    <- params[4]
	c       <- params[5]
	
	if (all(is.na(noise.targets))) { 
		noise.targets <- rep(sigma.n, length(noise.targets))
	}
	
	n <- length(y)
	n2 <- length(x.targets)
	y <- y - mean(y)
	
	# shortcuts
	bexp2 <- function(x) { .bexp(x, l, lmin, c) }
	s1 <- function(t1,t2) {    return (t1/bexp2(t1) - t2/bexp2(t2))  }
	s2 <- function(t1,t2) {    return ((t2-t2*exp(-c*t2))/bexp2(t2)^2 - (t1-t1*exp(-c*t1))/bexp2(t1)^2)  }
	s3 <- function(t1,t2) {    return ((t2*exp(-c*t2))/bexp2(t2)^2 - (t1*exp(-c*t1))/bexp2(t1)^2)  }
	s4 <- function(t1,t2) {    return ((l-lmin) * ( (t2^2*exp(-c*t2))/bexp2(t2)^2 - (t1^2*exp(-c*t1))/bexp2(t1)^2))  }
	
	# compute kernels
	K = .nsRBF.blocks(x, x.targets, sigma.f, l, lmin, c)
	K.tt = K$tt
	K.ts = K$ts
	K.st = K$st
	K.ss = K$ss
	#  K.tt    <- .nsRBF(x, x, sigma.f, l, lmin, c) # data cov
	#  K.ts    <- .nsRBF(x, x.targets, sigma.f, l, lmin, c) # cross
	#  K.st    <- .nsRBF(x.targets, x, sigma.f, l, lmin, c) # cross
	#  K.ss    <- .nsRBF(x.targets, x.targets, sigma.f, l, lmin, c) # background
	K.y     <- K.tt + diag(noise.targets[match(x, x.targets)])^2 # add stdev noise term
	K.y.inv <- solve(K.y)
	f.mean  <- K.st%*%K.y.inv%*%y
	f.cov   <- K.ss - K.st%*%K.y.inv%*%K.ts
	K.m     <- K.ss + diag(noise.targets)^2
	K.m.inv <- solve(K.m)
	
	# more stuff
	alphas <- K.m.inv %*% f.mean
	betas <- K.y.inv %*% K.ts
	
	# some other matrices
	B <- ((2*t(betas))%*%K.ts - t(betas)%*%K.tt%*%betas + 2*K.ss)
	S1 <- mat.or.vec(n2,n2)
	S2 <- mat.or.vec(n2,n2)
	S3 <- mat.or.vec(n2,n2)
	S4 <- mat.or.vec(n2,n2)
	for (i in 1:n2) {
		for (j in 1:n2) {
			t1 <- x.targets[i]
			t2 <- x.targets[j]
			S1[i,j] <- s1(t1,t2)
			S2[i,j] <- s2(t1,t2)
			S3[i,j] <- s3(t1,t2)
			S4[i,j] <- s4(t1,t2)
		}
	}
	
	# gradients
	gKf    <-  0.5 * sigma.f * B
	gKl    <- -2 * S1 * S2 * B
	gKlmin <- -2 * S1 * S3 * B
	gKc    <- -2 * S1 * S4 * B
	
	gf    <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% gKf))
	gn    <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% (2*sigma.n*diag(n2)))) + (alpha-1)/sigma.n - beta
	gl    <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% gKl))
	glmin <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% gKlmin))
	gc    <- 0.5* sum(diag((alphas %*% t(alphas) - K.m.inv) %*% gKc))
	
	print (c(gf,gn,gl,glmin,gc))
	
	return (c(gf,gn,gl,glmin,gc))
}


.randparams = function(defaults, lbounds, ubounds, n=5) {

	initials <- mat.or.vec(n,5)
	initials[1,] <- defaults
	for (i in 2:n) {
		initials[i,1] <- runif(1, lbounds[1], ubounds[1]) # f
		initials[i,2] <- runif(1, lbounds[2], ubounds[2]) # n
		initials[i,3] <- runif(1, lbounds[3], ubounds[3]) # l
		initials[i,4] <- runif(1, lbounds[4], ubounds[4]) # lmin
		initials[i,5] <- runif(1, lbounds[5], ubounds[5]) # c
		
		if (initials[i,4] > initials[i,3]) # swap
		{
			initials[i,3:4] = rev(initials[i,3:4])
		}
	}
	
	return(initials)
}


# expected optimization:
# use L-BFGS-B gradient ascent with analytical gradients
# do 
# .. 1 optimization with good initial values
# .. 4 with random initial values
# for a total of 5 rounds (difficult to optimize)
.optimize_expected <- function(x, y, x.targets, noise.targets, a, b, defaultparams, lbounds, ubounds, restarts=3) {
	initials = .randparams(defaultparams, lbounds, ubounds, restarts)
	
	params = defaultparams
	bestval <- -Inf
	for (i in 1:restarts) {
		res <- try(optim(par=initials[i,],
							  fn= .optim_exp,
							  gr= .optim_exp_gradients, lower=lbounds, upper=ubounds, method='L-BFGS-B', 
							  control=list(fnscale=-1), x=x, y=y, x.targets, noise.targets=noise.targets, alpha=a, beta=b))
		if (class(res) != 'try-error' && res$value > bestval) {
			params <- c(res$par[1], 1, res$par[2:4])
			bestval <- res$value
		}
	}
	
	return (params)
}

# dynamic optimization:
# use L-BFGS-B gradient ascent with analytical gradients
# do
# .. 5 optimizations with random initial values
# .. 1 with expected 'good' values
# .. 1 with first optimizing only likelihood, then 1 also model complexity
# for a total of 8 rounds
.optimize_dynamic <- function(x, y, x.targets, noise.targets, a, b, defaultparams, lbounds, ubounds, restarts=3) {
	initials = .randparams(defaultparams, lbounds, ubounds, restarts)


	params = defaultparams
	bestval <- -Inf
	for (i in 1:restarts) {
		res <- try(optim(par=initials[i,],
							  fn= .optim_dyn,
							  gr= .optim_dyn_gradients, lower=lbounds, upper=ubounds, method='L-BFGS-B', 
							  control=list(fnscale=-1,parscale=rep(0.01,5)), x=x, y=y, x.targets=x.targets, noise.targets=noise.targets, alpha=a, beta=b))
		if (class(res) != 'try-error' && res$value > bestval) {
			params <- res$par
			bestval <- res$value
		}
	}
	
	lambda = 0.2
	res <- try(optim(par=defaultparams,
						  fn= .optim_fit_dyn,
						  gr= .optim_fit_dyn_gradients, lower=lbounds, upper=ubounds, method='L-BFGS-B', 
						  control=list(fnscale=-1,parscale=rep(0.01,5)), x=x, y=y, x.targets=x.targets, noise.targets=noise.targets, alpha=a, beta=b, lambda=lambda))
	
	if (class(res) == 'try-error') {
		return (params)
	}
	
	res2 <- try(optim(par=res$par,
							fn= .optim_dyn,
							gr= .optim_dyn_gradients, lower=lbounds, upper=ubounds, method='L-BFGS-B', 
							control=list(fnscale=-1,parscale=rep(0.01,5)), x=x, y=y, x.targets=x.targets, noise.targets=noise.targets, alpha=a, beta=b))
	if (class(res) == 'try-error') {
		return (params)
	}
	
	if (res2$value > bestval) {
		params <- res2$par
	} 
	
	return(params)
}

# optimize the static model, easy
.optimize_static <- function(x, y, x.targets, noise.targets, a, b, defaultparams, lbounds, ubounds, restarts) {
	initials = .randparams(defaultparams, lbounds, ubounds, restarts)
	
	params = defaultparams
	bestval = -Inf
	for (i in 1:restarts) {
		res <- try(optim(par=initials[i,1:3],
							  fn= .optim_stat,
							  gr= .optim_stat_gradients,
							  lower=lbounds[1:3],
							  upper=ubounds[1:3],
							  method="L-BFGS-B", 
							  control=list(fnscale=-1), x=x,y=y, x.targets=x.targets, noise.targets=noise.targets, alpha=a, beta=b))
		if (class(res) != 'try-error' && res$value > bestval) {
			params <- res$par
			bestval <- res$value
		}
	}
	
	return (params)
}


