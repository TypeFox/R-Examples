# auxiliary functions for GPR
#


# normalise kernel matrix
.normkernel = function(K) {
	return(K / sqrt(diag(K)%*%t(diag(K))))
}

# gaussian kernel function
.gaussiankernel <- function(x1, x2, sf=1, l=1) {
	return(sf^2*exp(-(x1-x2)^2/l^2))
}

# bounded exponential function
.bexp <- function(x, max, min=0, curve=Inf) {
	return (max - (max-min)*exp(-curve*(x+0.0001)))
}

# triangular kernel matrix
.triangular <- function(X1, X2, dist=5) {
	K <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
	for (i in 1:nrow(K)) {
		for (j in 1:ncol(K)) {
			d = X1[i] - X2[j]
			K[i,j] <- max(c(0, (1 - abs(d)/dist)))
		}
	}
	
	return(K)
}

# linear kernel matrix
.uniform <- function(X1, X2, dist=5) {
	K <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
	for (i in 1:nrow(K)) {
		for (j in 1:ncol(K)) {
			if (abs(X1[i] - X2[j]) < dist) {
				K[i,j] = 1/(dist*2);
			}
		}
	}
	
	return(K)
}

# gaussian kernel matrix
.RBF <- function(X1, X2, sf=1, l=1) {
	K <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
	for (i in 1:nrow(K)) {
		for (j in 1:ncol(K)) {
			K[i,j] <- sf^2 * exp(-(X1[i]-X2[j])^2/l^2)
		}
	}
	
	return(K)
}

# derivative of gaussian kernel (matrix)
.RBF.deriv = function(X1,X2, sf=1, l=1) {
	Kd <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
	for (i in 1:nrow(Kd)) {
		for (j in 1:ncol(Kd)) {
			Kd[i,j] <- -sf^2 * 2*(X1[i]-X2[j])/l^2 * exp(-(X1[i]-X2[j])^2/l^2)
		}
	}
	
	return (Kd)
}

# second derivative of gaussian kernel (matrix)
.RBF.deriv2 = function(X1,X2, sf=1, l=1) {
	Kd <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
	for (i in 1:nrow(Kd)) {
		for (j in 1:ncol(Kd)) {
			Kd[i,j] <- sf^2 * 4*(X1[i]-X2[j])^2/l^4 * exp(-(X1[i]-X2[j])^2/l^2)
		}
	}
	
	return (Kd)
}

# non-stationary gaussian kernel matrix
.nsRBF <- function(X1, X2, sf=1, l=1,lmin=0,c=Inf) {
	K <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
	for (i in 1:nrow(K)) {
		for (j in 1:ncol(K)) {
			K[i,j] <- sf^2 * exp(-( abs(X1[i])/.bexp(X1[i],l,lmin,c) - abs(X2[j])/.bexp(X2[j],l,lmin,c) )^2)
		}
	}
	
	return(K)
}

# derivative of non-stationary gaussian kernel (matrix)
.nsRBF.deriv <- function(X1, X2, sf=1, l=1,lmin=0,c=Inf) {
	Kd <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
	for (i in 1:nrow(Kd)) {
		for (j in 1:ncol(Kd)) {
			Kd[i,j] <- -sf^2 * 2*( X1[i]/.bexp(X1[i],l,lmin,c) - X2[j]/.bexp(X2[j],l,lmin,c) )/.bexp(X1[i],l,lmin,c) * exp(-( abs(X1[i])/.bexp(X1[i],l,lmin,c) - abs(X2[j])/.bexp(X2[j],l,lmin,c) )^2)
		}
	}
	
	return(Kd)
}

# second derivative of non-stationary gaussian kernel (matrix)
.nsRBF.deriv2 <- function(X1, X2, sf=1, l=1,lmin=0,c=Inf) {
	Kd <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
	for (i in 1:nrow(Kd)) {
		for (j in 1:ncol(Kd)) {
			Kd[i,j] <- sf^2 * 4*( X1[i]/.bexp(X1[i],l,lmin,c) - X2[j]/.bexp(X2[j],l,lmin,c) )^2/.bexp(X1[i],l,lmin,c)^2 * exp(-( abs(X1[i])/.bexp(X1[i],l,lmin,c) - abs(X2[j])/.bexp(X2[j],l,lmin,c) )^2)
		}
	}
	
	return(Kd)
}

# gaussian kernel matrix, decomposed into four blocks
.RBF.blocks = function(x1, x2, sigma.f, l) {
	K = .RBF(c(x1,x2), c(x1,x2), sigma.f, l)
	
	n1 = length(x1)
	n2 = length(x1) + length(x2)
	
	k.tt = K[1:n1,      1:n1]
	k.ts = K[1:n1,      (n1+1):n2]
	k.st = K[(n1+1):n2, 1:n1]
	k.ss = K[(n1+1):n2, (n1+1):n2]
	
	return(list('tt'=k.tt, 'ts'=k.ts, 'st'=k.st, 'ss'=k.ss))
}

# non-stationary gaussian kernel matrix, decomposed into four blocks
.nsRBF.blocks = function(x1, x2, sigma.f, l, lmin, c) {
	K = .nsRBF(c(x1,x2), c(x1,x2), sigma.f, l, lmin, c)
	
	n1 = length(x1)
	n2 = length(x1) + length(x2)
	
	k.tt = K[1:n1,      1:n1]
	k.ts = K[1:n1,      (n1+1):n2]
	k.st = K[(n1+1):n2, 1:n1]
	k.ss = K[(n1+1):n2, (n1+1):n2]
	
	return(list('tt'=k.tt, 'ts'=k.ts, 'st'=k.st, 'ss'=k.ss))
}

.learnnoise = function(x, y, x.targets) {
	vars = tapply(y, x, FUN=var) # compute variances
	vars = vars[!is.na(vars)] # remove NA's
	noise = .interpolate(names(vars), vars, x.targets)  # interpolate noises
	noise = sqrt(noise)
	
	return(noise)
}

.interpolate = function(x, vals, targets) {
	spar = 0.35
	good = !is.na(vals) # remove NA's
	x = x[good]
	vals = vals[good]
	
	# not enough replicates
	if (length(vals) < 2) {
		#    print('Not enough replicates to learn observational noise\nGive noise with noise.ctrl and noise.case')
		noise = rep(NA, length(targets))
		return (noise)
	}
	
	noise = predict(smooth.spline(x, vals, spar=spar), targets)$y
	noise[noise < min(vals)] = min(vals) # don't let the curve go under
	
	return(noise)
}

.handlenoise = function(x, y, x.targets, noise, nsnoise) {
	n = length(x)
	n.targets = length(x.targets)
	
	# we are given noises for observations -> interpolate them
	if (nsnoise==TRUE && length(noise) == n) {
		noise = .interpolate(x, noise, x.targets)
	}
	
	# we are not given noises -> learn from replicates
	if (nsnoise==TRUE && (is.null(noise) || is.na(noise))) {
		noise = .learnnoise(x, y, x.targets)
	}
	
	if (is.null(noise) || is.na(noise)) {  # we are not given noise -> set as NA and learn global noise level later
		noise.targets = rep(NA, n.targets)
		noise.obs = rep(NA, n)
	} else if (length(noise) == 1) {   # we are given static noise -> use for all times
		noise.targets = rep(noise, n.targets)
		noise.obs = rep(noise, n)
	} else if (length(noise) == n) {   # we are given noise.obs -> use same for targets
		noise.targets = rep(NA, n.targets)
		noise.targets[match(x, x.targets)] = noise
		noise.obs = noise
	} else if (length(noise) == n.targets) {   # we are given noise.targets -> use same for obs
		noise.targets = noise
		noise.obs = noise.targets[match(x, x.targets)]
	}
	
	return(list(noise.obs=noise.obs,noise.targets=noise.targets))
}




