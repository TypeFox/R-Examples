estim.gsm <- function(y, J, G = 100, M = 600, a, b, alpha, init = list(rep(1 / J, J), NA, rep(1, N))) {
	N <- length(y)
	y.grid <- seq(min(y)*.66, max(y)*1.5, length = G)
	
	lbl <- matrix(NA, nrow = N, ncol = M)
	wgt <- matrix(NA, nrow = M, ncol = J)
	
	w <- init[[1]]
	theta <- rep(init[[2]], M)
	x <- init[[3]]
	logfdens <- matrix(NA, J, G)
	fdens <- matrix(NA, M, G)
	
	for (m in 1:M) {
		for (i in 1:N) {
			temp <- log(seq(a + sum(x) - x[i], a + sum(x) - x[i] + J - 1))
			pi <- (1:J - 1)*log(y[i]) - lgamma(1:J) + cumsum(temp) - (1:J)*log(b + sum(y))
			pi <- (w*exp(pi)) / sum(w*exp(pi))
			x[i] <- sample(1:J, 1, prob = pi)
			lbl[i, m] <- x[i]
		}
		
		x.counts <- table(factor(x, levels = 1:J))
		w <- rdirichlet(1, x.counts + alpha)
		theta[m] <- rgamma(1, sum(x) + a, rate = sum(y) + b)
		for (j in 1:J) {
			logfdens[j,] <- log(w[j]) + (j - 1)*log(y.grid) - theta[m]*y.grid - lgamma(j) + j*log(theta[m])
		}
		fdens[m, ] <- apply(exp(logfdens), 2, sum)
		wgt[m, ] <- w
		if (m / 100 == round(m / 100)) print(paste("simulation", m, "/", M))
	}

	out <- new("gsm", fdens = fdens, theta = theta, weight = wgt, label = lbl, data = y)
	return(out)
}

estim.gsm_theta <- function(y, J, G = 100, M = 600, a, b, alpha, init = list(rep(1 / J, J), J / max(y), rep(1, N))) {
	N <- length(y)
	y.grid <- seq(min(y)*.66, max(y)*1.5, length = G)
	
	lbl <- matrix(NA, nrow = N, ncol = M)
	wgt <- matrix(NA, nrow = M, ncol = J)
	
	w <- init[[1]]
	theta <- rep(init[[2]], M + 1)
	x <- init[[3]]
	logfdens <- matrix(NA, J, G)
	fdens <- matrix(NA, M, G)
	
	for (m in 1:M) {
		for (i in 1:N) {
			pi <- (1:J - 1)*log(y[i]) - theta[m]*y[i] - lgamma(1:J) + (1:J)*log(theta[m])
			pi <- (w*exp(pi)) / sum(w*exp(pi))
			x[i] <- sample(1:J, 1, prob = pi)
			lbl[i, m] <- x[i]
		}
		
		x.counts <- table(factor(x, levels=1:J))
		w <- rdirichlet(1, x.counts + alpha)
		theta[m + 1] <- rgamma(1, sum(x) + a, rate = sum(y) + b)
		for (j in 1:J) {
			logfdens[j,] <- log(w[j]) + (j - 1)*log(y.grid) - theta[m + 1]*y.grid - lgamma(j) + j*log(theta[m + 1])
		}
		fdens[m, ] <- apply(exp(logfdens), 2, sum)
		wgt[m, ] <- w
		if (m / 100 == round(m / 100)) print(paste("simulation", m, "/", M))
	}

	out <- new("gsm", fdens = fdens, theta = theta, weight = wgt, label = lbl, data = y)
	return(out)
}

allcurves.q <- function(post, perc) {
	n <- dim(post)[2]
	temp <- numeric(n)
	for (i in 1:n) temp[i] <- quantile(post[, i], perc) # to be replaced with 'apply(post, 2, quantile, probs = perc)'
	return(temp)
}

dgsm <- function(x, weight, rateparam) {
	numcomp <- length(weight)
	mixcomp <- matrix(NA, nrow = length(x), ncol = numcomp)
	for (i in 1:numcomp) mixcomp[ , i] <- dgamma(x, shape = i, rate = rateparam)
	dens <- as.numeric(mixcomp%*%weight)
	return(dens)
}

pgsm <- function (q, weight, rateparam, lower.t = TRUE) {
	numcomp <- length(weight)
	mixcomp <- matrix(NA, nrow = length(q), ncol = numcomp)
	for (i in 1:numcomp) mixcomp[, i] <- pgamma(q, shape = i, rate = rateparam, lower.tail = lower.t)
	cumprob <- as.numeric(mixcomp %*% weight)
	return(cumprob)
}

rgsm <- function(n, weight, rateparam) {
	J <- length(weight)
	rmixg <- vector(length = n)
	tmp.rmixg <- matrix(NA, nrow = J, ncol = n)
	tmp.lbl <- rmultinom(n, 1, weight)
	for (i in 1:J) {
		tmp.rmixg[i,] <- rgamma(n, shape = i, rate = rateparam)
		rmixg <- rmixg + tmp.lbl[i,]*tmp.rmixg[i,]
	}
	return(rmixg)
}

qgsm <- function(p, x = NULL, weight, rateparam, alpha = .05, br = c(0, 1000), lower.t = TRUE) {
	if (any(p < 0 | p > 1)) stop("probabilities outside the unit interval.")
	if (alpha > .5) stop("alpha has to be strictly smaller than .5.")
	
	F_inv_tmp <- function(p, weight, rateparam, br, lower.t) {
		G <- function(x) {
			pgsm(x, weight, rateparam, lower.t) - p
		}
		return(uniroot(f = G, interval = br)$root)
	}
	
	qgsm_tmp <- Vectorize(F_inv_tmp, "p")
	
	n <- length(p)
	out <- numeric(n)
	p.indx <- which((p >= alpha) & (p <= (1 - alpha)))
	if (!is.null(x)) {
		if (length(p[p.indx]) > 0) {
			p.tmp <- pgsm(q = x, weight = weight, rateparam = rateparam, lower.t = lower.t)
			cdf.gsm <- sortedXyData(x, p.tmp)
			q1 <- sapply(p[p.indx], FUN = NLSstClosestX, xy = cdf.gsm, simplify = TRUE)
			if (length(p[-p.indx]) > 0) {
				q2 <- qgsm_tmp(p = p[-p.indx], weight = weight, rateparam = rateparam, br = br, lower.t = lower.t)
				out[p.indx] <- q1
				out[-p.indx] <- q2 }
			else {
				out <- q1 } }
		else {
			out <- qgsm_tmp(p = p, weight = weight, rateparam = rateparam, br = br, lower.t = lower.t) } }
	else {
		out <- qgsm_tmp(p = p, weight = weight, rateparam = rateparam, br = br, lower.t = lower.t) }
	return(out)
}
