llDiffT <- function(dat, thetai, cand, deltaj, gamma, sigma2) {
	y <- dat$y
	N <- dat$N
	groupi <- dat$groupi
	groupj <- dat$groupj
	yi <- tapply(y, groupi, sum)
	logdiffs <- log(1 + exp(cand[groupi] + deltaj[groupj])) - log(1 + exp(thetai[groupi] + deltaj[groupj]))
	out <- as.numeric(yi*(cand - thetai) - tapply(N*logdiffs, groupi, sum) - (cand^2 - thetai^2)/(2*sigma2))
	return(out)
}

llDiffD <- function(dat, deltaj, cand, thetai, gamma, tau2) {
	y <- dat$y
	N <- dat$N
	groupi <- dat$groupi
	groupj <- dat$groupj
	y.j <- tapply(y, groupj, sum)
	logdiffs <- log(1 + exp(thetai[groupi] + cand[groupj])) - log(1 + exp(thetai[groupi] + deltaj[groupj]))
	out <- as.numeric(y.j*(cand - deltaj) - tapply(N*logdiffs, groupj, sum) - ((cand - gamma)^2 - (deltaj - gamma)^2)/(2*tau2))
	return(out)
}

bhm.mcmc <- function(dat, nsim = 2000, burnin = 500, scale.factor = 1, adaptive.int = 100, adaptive.max = 1000,
		prior = NULL, init = NULL, tuneD = NULL, tuneT = NULL) {
	
	colVs <- function(x) {
		out <- (colSums(x*x)-nrow(x)*colMeans(x)^2)/(nrow(x)-1)
		return(out)
	}
	
	if (class(dat) != "mederrData") {
		stop("dat argument not of class 'mederrData'. Type '?mederrData' for help.")
	}
	
	data <- dat@data
	y <- data$y
	N <- data$N
	groupi <- factor(data$groupi)
	groupj <- factor(data$groupj)
	
	ni <- nlevels(groupi)
	nj <- nlevels(groupj)
	M <- burnin + nsim
	
	if (is.null(tuneT)) tuneT <- rep(.5, ni)
	if (is.null(tuneD)) tuneD <- rep(.4, nj)
	
	if (is.null(init)) {
		tab1 <- rowsum(data[, 1:2], groupi)
		ti <- log(tab1$y/(tab1$N - tab1$y))
		ti[abs(ti) > 1000] <- mean(ti[abs(ti) < 1000])
		tab2 <- rowsum(data[, 1:2], groupj)
		dj <- log(tab2$y/(tab2$N - tab2$y)) - mean(ti)
		dj[abs(dj) > 1000] <- 0
		init <- list(thetai = ti, deltaj = dj, gamma = mean(ti), sigma2 = var(ti), tau2 = var(dj))
	}
	thetai <- init$thetai
	deltaj <- init$deltaj
	gamma <- init$gamma
	sigma2 <- init$sigma2
	tau2 <- init$tau2
	
	if (is.null(prior)) { 
		prior <- list(g = -4, G = 1000, a1 = .001, b1 = .001, a2 = .001, b2 = .001)
	}
	
	thetai.chain <- matrix(nrow = M, ncol = ni)
	deltaj.chain <- matrix(nrow = M, ncol = nj)
	gamma.chain <- sigma2.chain <- tau2.chain <- numeric(length = M)
	
	counti <- rep(0, ni)
	countj <- rep(0, nj)
	
	for(s in 1:M) {
		tau2 <- 1/rgamma(1, prior$a2 + nj/2, prior$b2 + .5*sum((deltaj-gamma)^2))
		sigma2 <- 1/rgamma(1, prior$a1 + ni/2, prior$b1 + .5*sum(thetai^2))
		V <- 1/(1/prior$G + nj/tau2)
		gamma <- rnorm(1, mean = V * (prior$g/prior$G + sum(deltaj)/tau2), sd = sqrt(V))
		
		if ((s%%adaptive.int == 0) & (s <= adaptive.max)) {
			tuneD <- sqrt(scale.factor*colVs(deltaj.chain[1:(s-1), ]))
			tuneT <- sqrt(scale.factor*colVs(thetai.chain[1:(s-1), ]))
		}

		cand <- rnorm(n = nj, mean = deltaj, sd = tuneD)
		accept <- llDiffD(data, deltaj, cand, thetai, gamma, tau2) > log(runif(1))
		if (s > burnin) countj <- countj + accept
		deltaj <- ifelse(accept, cand, deltaj)

		cand <- rnorm(n = ni, mean = thetai, sd = tuneT)
		accept <- llDiffT(data, thetai, cand, deltaj, gamma, sigma2) > log(runif(1))
		if (s > burnin) counti <- counti + accept
		thetai <- ifelse(accept, cand, thetai)
		
		thetai.chain[s, ] <- thetai
		deltaj.chain[s, ] <- deltaj
		gamma.chain[s] <- gamma
		sigma2.chain[s] <- sigma2
		tau2.chain[s] <- tau2

		if ((s/100 == round(s/100))) print(paste("iter ", s, " / ", M, sep = ""), quote = FALSE)
	}
	out <- new("mederrFit", thetai = thetai.chain[(burnin+1):M, ], deltaj = deltaj.chain[(burnin+1):M, ], gamma = gamma.chain[(burnin+1):M],
			sigma2 = sigma2.chain[(burnin+1):M], tau2 = tau2.chain[(burnin+1):M], p.acc.i = counti/nsim, p.acc.j = countj/nsim, tune.theta = tuneT,
			tune.delta = tuneD, k = Inf, eta = 1)
	return(out)
}

dst <- function(x, sigma, k, eta) {
	out <- (2/(eta + 1/eta)/sigma)*ifelse(x >= 0, dt(x/eta/sigma, df = k), dt(x*eta/sigma, df = k))
	return(out)
}

logp <- function(theta, deltaj, sigma2, i, k, eta, dat) {

	if (class(dat) != "mederrData") {
		stop("dat argument not of class 'mederrData'. Type '?mederrData' for help.")
	}
	
	if (is.vector(deltaj)) deltaj <- t(deltaj)
	p <- theta + deltaj
	out <- - colSums(dat@data$y[dat@data$groupi == i]*t(p) - dat@data$N[dat@data$groupi == i]*log(1+exp(t(p)))) - log(dst(theta, sqrt(sigma2), k, eta))
	return(out)
}

bhm.resample <- function(model, dat, p.resample = 0.1, k, eta) {

	if (class(model) != "mederrFit") {
		stop("Model argument not of class 'mederrFit'. Type '?mederrFit' for help.")
	}
	if (class(dat) != "mederrData") {
		stop("dat argument not of class 'mederrData'. Type '?mederrData' for help.")
	}

	deltaj.hat <- rowMeans(model@deltaj)
	sigma2.hat <- mean(model@sigma2)
	
	grd <- list("k" = k, "eta" = eta)
	nk <- length(grd$k)
	neta <- length(grd$eta)
	
	M <- nrow(model@thetai)
	L <- floor(M*p.resample)
	ni <- ncol(model@thetai)
	nj <- ncol(model@deltaj)
	
	log.ir <- array(dim = c(nk, neta, M))
	samp <- array(dim = c(nk, neta, L))
	A <- t.new <- array(dim = c(nk, neta, ni))
	p.old <- p.new <- matrix(nrow = M, ncol = ni)
	modes <- vector(mode = "list", length = ni)
	
	ow <- options("warn")
	options(warn = -1)
	for(i in 1:ni) {
		modes[[i]] <- nlm(logp, p = 0, hessian = TRUE, deltaj = deltaj.hat, sigma2 = sigma2.hat, i = i, k = Inf, eta = 1, dat)
		p.old[, i] <- logp(model@thetai[, i], deltaj = model@deltaj, sigma2 = model@sigma2, i = i, k = Inf, eta = 1, dat) 
	}
	t.old <- unlist(sapply(modes, function(x) x["estimate"]))
	names(t.old) <- NULL
	h.old <- unlist(sapply(modes, function(x) x["hessian"]))
	
	for(x in 1:nk) {
		k <- grd$k[x]
		
		for(y in 1:neta) {
			eta <- grd$eta[y]
			
			for(i in 1:ni) { 
				modes[[i]] <- nlm(logp, p = 0, hessian = TRUE, deltaj = deltaj.hat, sigma2 = sigma2.hat, i = i, k = k, eta = eta, dat)
			}
			t.new[x, y, ] <- unlist(sapply(modes, function(x) x["estimate"]))
			h.new <- unlist(sapply(modes, function(x) x["hessian"]))
			
			A[x, y, ] <- sqrt(h.old / h.new)
			lf <- t.new[x, y, ] + A[x, y, ]*(model@thetai - t.old)
			
			for(i in 1:ni) {
				p.new[, i] <- logp(lf[, i], deltaj = model@deltaj, sigma2 = model@sigma2, i = i, k = k, eta = eta, dat)
			}
			
			log.ir[x, y, ] <- sum(log(A[x, y, ])) + rowSums(p.old) - rowSums(p.new)
			lir.fix <- log.ir[x, y, ] + 500
			ir <- exp(lir.fix)
			samp[x, y, ] <- sample(M, L, prob = ir/sum(ir))
			
			print(paste("step ", (y + (x - 1)*neta), " / ", (nk*neta), sep = ""), quote = FALSE)
		}
	}
	options(ow)
	out <- new("mederrResample", log.ir = log.ir, samp = samp, A = A, t.new = t.new, t.old = t.old, grd = grd)
	return(out)
}

logunpost <- function(resample) {
	if (class(resample) != "mederrResample") {
		stop("Resample argument not of class 'mederrResample'. Type '?mederrResample' for help.")
	}
	
	out <- log(apply(exp(resample@log.ir), 1:2, sum))
	rownames(out) <- resample@grd$k
	colnames(out) <- resample@grd$eta

	return(out)
}

bhm.constr.resamp <- function(model, resample, k = Inf, eta = 1) {
	if (class(model) != "mederrFit") {
		stop("model argument not of class 'mederrFit'. Type '?mederrFit' for help.")
	}
	if (class(resample) != "mederrResample") {
		stop("Resample argument not of class 'mederrResample'. Type '?mederrResample' for help.")
	}
	
	i <- which(resample@grd$k == k)
	j <- which(resample@grd$eta == eta)
	ss <- resample@samp[i, j, ]
	
	thetai.resamp <- resample@t.new[i, j, ] + resample@A[i, j, ]*(model@thetai[ss, ] - resample@t.old)
	deltaj.resamp <- model@deltaj[ss, ]
	sigma2.resamp <- model@sigma2[ss]
	tau2.resamp <- model@tau2[ss]
	gamma.resamp <- model@gamma[ss]

	out <- new("mederrFit", thetai = thetai.resamp, deltaj = deltaj.resamp, gamma = gamma.resamp, sigma2 = sigma2.resamp, tau2 = tau2.resamp,
			p.acc.i = model@p.acc.i, p.acc.j = model@p.acc.j, tune.theta = model@tune.theta, tune.delta = model@tune.delta, k = k, eta = eta)
	return(out)
}

bayes.rank <- function(model) {
	if (class(model) != "mederrFit") {
		stop("Model argument not of class 'mederrFit'. Type '?mederrFit' for help.")
	}
	
	theta <- model@thetai
	ni <- ncol(theta)
	M <- nrow(theta)
	rank <- numeric(length = ni)
	for(i in 1:ni) {
		rank[i] <- sum((theta[,i] - theta) >= 0)/M
	}
	out <- rank
	return(out)
}

post.rep <- function(model, dat) {
	if (class(model) != "mederrFit") {
		stop("Model argument not of class 'mederrFit'. Type '?mederrFit' for help.")
	}
	if (class(dat) != "mederrData") {
		stop("dat argument not of class 'mederrData'. Type '?mederrData' for help.")
	}
	
	theta <- model@thetai
	delta <- model@deltaj
	y <- dat@data$y
	N <- dat@data$N
	groupi <- dat@data$groupi
	groupj <- dat@data$groupj
	
	ni <- ncol(theta)
	nj <- ncol(delta)
	M <- nrow(theta)
	
	yi.obs <- rowsum(y, groupi)
	yj.obs <- rowsum(y, groupj)
	
	yi.rep <- matrix(0, nrow = M, ncol = ni)
	yj.rep <- matrix(0, nrow = M, ncol = nj)
	for(s in 1:M) {
		l <- theta[s, ][groupi] + delta[s, ][groupj]
		y.rep <- rbinom(ni*nj, N, exp(l)/(1+exp(l)))
		yi.rep[s, ] <- rowsum(y.rep, groupi)
		yj.rep[s, ] <- rowsum(y.rep, groupj)
	}
	out <- list("yi.obs" = yi.obs, "yj.obs" = yj.obs, "yi.rep" = yi.rep, "yj.rep" = yj.rep)
	return(out)
}

p.value <- function(reps) {
	yi.obs <- reps$yi.obs
	yj.obs <- reps$yj.obs
	yi.rep <- reps$yi.rep
	yj.rep <- reps$yj.rep

	ni <- length(yi.obs)
	nj <- length(yj.obs)

	plessi <- pgreati <- numeric(length = ni)
	for(i in 1:ni) {
		plessi[i] <- mean(yi.rep[,i] <= yi.obs[i])
		pgreati[i] <- mean(yi.rep[,i] >= yi.obs[i])
	}
	plessj <- pgreatj <- numeric(length = nj)
	for(j in 1:nj) {
		plessj[j] <- mean(yj.rep[, j] <= yj.obs[j])
		pgreatj[j] <- mean(yj.rep[, j] >= yj.obs[j])
	}

	pvaluei <- pmin(plessi, pgreati)
	pvaluej <- pmin(plessj, pgreatj)
	out <- list("Yi" = pvaluei, "Yj" =  pvaluej)
	return(out)
}

dnegbinom <- function(x, alpha, prob, log.p = FALSE) {
	logdens <- lgamma(x + alpha) - lgamma(alpha) - lfactorial(x) + alpha*log(prob) + x*log(1 - prob)
	if (log.p) return(logdens) else return(exp(logdens))
}

negbinom.loglik <- function(theta, N, E) {
	loglik <- dnegbinom(N, theta[1], 1/(1 + E/theta[2]), log.p = TRUE)
	return(sum(loglik))
}

dmixnegbinom <- function(x, theta, E, log.p = FALSE) {
	logf1 <- dnegbinom(x, theta[1], 1/(1 + E/theta[2]), log.p = TRUE)
	logf2 <- dnegbinom(x, theta[3], 1/(1 + E/theta[4]), log.p = TRUE)
	logdens <- log(theta[5]*exp(logf1) + (1 - theta[5])*exp(logf2))
	if (log.p) return(logdens) else return(exp(logdens))
}

rmixnegbinom <- function (n, theta, E) {
	rmix <- vector(length = n)
	tmp.rmix <- matrix(NA, nrow = 2, ncol = n)
	tmp.lbl <- matrix(NA, nrow = 2, ncol = n)
	tmp.rmix[1, ] <- rnbinom(n, size = theta[1], prob = 1/(1 + E/theta[2]))
	tmp.rmix[2, ] <- rnbinom(n, size = theta[3], prob = 1/(1 + E/theta[4]))
	p <- c(theta[5], 1 - theta[5])
	tmp.lbl <- rmultinom(n, 1, p)
	for (i in 1:2) rmix <- rmix + tmp.lbl[i, ]*tmp.rmix[i, ]
	return(rmix)
}

mixnegbinom.loglik <- function(theta, N, E) {
	loglik <- dmixnegbinom(N, theta, E, log.p = TRUE)
	return(sum(loglik))
}

EBGM <- function(eb.result){
	
	N <- eb.result$N
	E <- eb.result$E
	theta <- eb.result$theta.hat
	prior <- eb.result$prior
	
	Q_n  <- function(theta, N, E){
		logf1 <- dnegbinom(N, theta[1], 1/(1 + E/theta[2]), log.p = TRUE)
		logf2 <- dnegbinom(N, theta[3], 1/(1 + E/theta[4]), log.p = TRUE)
		logdens <- log(theta[5]*exp(logf1) + (1 - theta[5])*exp(logf2))
		logQ_n <- log(theta[5]) + logf1 - logdens
		return(exp(logQ_n))
	}
	
	if (prior == "mixgamma") {
		a1 <- theta[1]
		b1 <- theta[2]
		a2 <- theta[3]
		b2 <- theta[4]
		P <- theta[5]
		Qn <- Q_n(theta, N, E)
		E_log.lambda <- Qn*(digamma(a1 + N) - log(b1 + E)) + (1 - Qn)*(digamma(a2 + N) - log(b2 + E))
	}
	else if (prior == "gamma") {
		a <- theta[1]
		b <- theta[2]
		E_log.lambda <- digamma(a + N) - log(b + E)
	}
	EBlog2 <- E_log.lambda/log(2)
	EBGM <- 2^EBlog2
	return(EBGM)
}

mixnegbinom.em <- function(dat, theta0, maxiter = 50000, toler = 0.01, se = TRUE, stratified = FALSE) {

	d.loglik.alpha <- function(alpha, N, E, theta) {
		n <- length(N)
		fnc <- log(1 + E/theta[n + 1]) - digamma(alpha + N) + digamma(alpha)
		return(sum(fnc*theta[1:n]))
	}

	d.loglik.beta <- function(beta, N, E, theta) {
		n <- length(N)
		fnc <- 1/(E + beta)*(N - theta[n + 1]*E/beta)
		return(sum(fnc*theta[1:n]))
	}
	
	if (class(dat) != "mederrData") {
		stop("dat argument not of class 'mederrData'. Type '?mederrData' for help.")
	}
	
	if (dat@size > 100000) cat("Preparing the data...")
	data <- dat@data
	err.dim <- dat@numi
	harm.dim <- 2
	hosp.dim <- dat@numj

	dat.sort <- data[sort(data$groupj, index.return = TRUE)$ix, ]
	tbl <- array(NA, c(harm.dim, err.dim, hosp.dim))
	for (j in 1:hosp.dim) tbl[, , j] <- t(dat.sort[(1 + err.dim*(j - 1)):(err.dim*j), 1:2])
	tbl[2, , ] <- tbl[2, , ] - tbl[1, , ]

	marg.err <- apply(tbl, c(1, 3), sum)
	marg.harm <- apply(tbl, c(2, 3), sum)
	joint.err.harm <- apply(tbl, c(1, 2), sum)
	marg.errharm <- colSums(marg.harm)

	if (stratified) {
		ind.joint.err.harm <- matrix(as.vector(t(marg.err))/rep(marg.errharm, harm.dim), nrow = harm.dim, ncol = hosp.dim, byrow = TRUE) %*% t(marg.harm)
	}
	else {
		marg.err <- colSums(joint.err.harm)
		marg.harm <- rowSums(joint.err.harm)
		total <- sum(joint.err.harm)
		ind.joint.err.harm <- marg.harm%*%t(marg.err)/total
	}

	x <- as.vector(t(joint.err.harm))
	E <- as.vector(t(ind.joint.err.harm))
	if (dat@size > 100000) cat("...done!\n")

	n <- length(x)
	error <- 1
	j <- 0
	alpha1 <- theta0[1]
	beta1 <- theta0[2]
	alpha2 <- theta0[3]
	beta2 <- theta0[4]
	pi.j <- c(theta0[5], 1 - theta0[5])
	theta.eb <- theta0
	
	dnbs <- rbind(t(dnegbinom(x, alpha1, 1/(1 + E/beta1))), t(dnegbinom(x, alpha2, 1/(1 + E/beta2))))
	pi.ij0 <- diag(pi.j) %*% dnbs
	pi.ij <- matrix(NA, nrow = 2, ncol = n)
	for (i in 1:n) pi.ij[,i] <- pi.ij0[, i]/colSums(pi.ij0)[i]
	
	while((error > toler) & (j < maxiter)) {
		oldp <- pi.ij
		olda1 <- alpha1
		olda2 <- alpha2
		oldb1 <- beta1
		oldb2 <- beta2
		oldll <- mixnegbinom.loglik(c(olda1, oldb1, olda2, oldb2, pi.j[1]), x, E)
		oldscore <- mixnegbinom.score(c(alpha1, beta1, alpha2, beta2, pi.j[1]), x, E)
		if (j == 0) print(paste("Iteration 1: ||score|| = ", formatC(sqrt(sum(oldscore^2)), digits = 4, format = "f"), " - loglik = ",
							formatC(oldll, digits = 4, format = "f"), sep = ""), quote = FALSE)
		
		theta <- c(oldp[1, ], oldb1)
		alpha1 <- BBsolve(olda1, d.loglik.alpha, N = x, E = E, theta = theta, quiet = TRUE)$par
		theta <- c(oldp[1, ], alpha1)
		beta1 <- BBsolve(oldb1, d.loglik.beta, N = x, E = E, theta = theta, quiet = TRUE)$par
		theta <- c(oldp[2, ], oldb2)
		alpha2 <- BBsolve(olda2, d.loglik.alpha, N = x, E = E, theta = theta, quiet = TRUE)$par
		theta <- c(oldp[2, ], alpha2)
		beta2 <- BBsolve(oldb2, d.loglik.beta, N = x, E = E, theta = theta, quiet = TRUE)$par
		pi.j <- apply(oldp, 1, sum)/n

		theta.eb <- c(alpha1, beta1, alpha2, beta2, pi.j[1])
		if (any(theta.eb <= 0) | theta.eb[5] > 1) {
			stop("Initial values not valid. Try a different set. Type '?mixnegbinom.em' for help.")
		}

		loglik <- mixnegbinom.loglik(theta.eb, x, E)
		score <- mixnegbinom.score(theta.eb, x, E)
		
		dnbs <- rbind(t(dnegbinom(x, alpha1, 1/(1 + E/beta1))), t(dnegbinom(x, alpha2, 1/(1 + E/beta2))))
		pi.ij0 <- diag(pi.j) %*% dnbs
		pi.ij <- matrix(NA, nrow = 2, ncol = n)
		for (i in 1:n) pi.ij[, i] <- pi.ij0[, i]/colSums(pi.ij0)[i]
		
		error <- abs(oldll - loglik)
		j <- j + 1
		if ((j/100 == round(j/100))) print(paste("iter ", j, " / ", maxiter, " - error = ", formatC(error, digits = 5, format = "f"), " > ",
							toler, " - ||score|| = ", formatC(sqrt(sum(score^2)), digits = 4, format = "f"), " - loglik = ",
							formatC(loglik, digits = 4, format = "f"), sep = ""), quote = FALSE)
	}
	print(paste("Iteration ", (j - 1)," (last): ||score|| = ", formatC(sqrt(sum(score^2)), digits = 4, format = "f"), " - loglik = ",
					formatC(loglik, digits = 4, format = "f"), sep = ""), quote = FALSE)

	if (se) {
		hess.lower <- genD(mixnegbinom.loglik, theta.eb, N = x, E = E)$D[, 6:20]
		hess <- matrix(NA, nrow = 5, ncol = 5)
		for (i in 1:5) {
			hess[i, (1:i)] <- hess.lower[(1 + sum(0:(i - 1))):sum(1:i)]
			hess[(1:i), i] <- hess[i, (1:i)]
		}
		se <- sqrt(diag(solve(-hess)))
		out <- list(theta.hat = theta.eb, final.err = error, final.ll = loglik, final.score = score, num.iter = (j - 1), se = se, N = x, E = E,
				prior = "mixgamma")
	}
	else {
		out <- list(theta.hat = theta.eb, final.err = error, final.ll = loglik, final.score = score, num.iter = (j - 1), se = NA, N = x, E = E,
				prior = "mixgamma")
	}

	return(out)
}

mixnegbinom.score <- function(theta, N, E) {
	n <- length(N)
	den <- dmixnegbinom(N, theta, E)
	fnc <- matrix(NA, nrow = n, ncol = 5)
	score <- rep(NA,5)
	fnc[, 1] <- theta[5]*dnegbinom(N, theta[1], 1/(1 + E/theta[2]))*(-log(1 + E/theta[2]) + digamma(theta[1] + N) - digamma(theta[1]))
	fnc[, 2] <- theta[5]*dnegbinom(N, theta[1], 1/(1 + E/theta[2]))*(-1/(E + theta[2])*(N - theta[1]*E/theta[2]))
	fnc[, 3] <- (1 - theta[5])*dnegbinom(N, theta[3], 1/(1 + E/theta[4]))*(-log(1 + E/theta[4]) + digamma(theta[3] + N) - digamma(theta[3]))
	fnc[, 4] <- (1 - theta[5])*dnegbinom(N, theta[3], 1/(1 + E/theta[4]))*(-1/(E + theta[4])*(N - theta[3]*E/theta[4]))
	fnc[, 5] <- dnegbinom(N, theta[1], 1/(1 + E/theta[2])) - dnegbinom(N, theta[3], 1/(1 + E/theta[4]))
	score <- apply(fnc/den, 2, sum)
	return(score)
}

negbinom.em <- function(dat, theta0, maxiter = 50000, toler = 0.01, se = TRUE, stratified = FALSE) {

	d.loglik.alpha <- function(alpha, N, E, theta) {
		fnc <- log(1 + E/theta) - digamma(alpha + N) + digamma(alpha)
		return(sum(fnc))
	}	
	d.loglik.beta <- function(beta, N, E, theta) {
		fnc <- 1/(E + beta)*(N - theta*E/beta)
		return(sum(fnc))
	}
	
	if (class(dat) != "mederrData") {
		stop("dat argument not of class 'mederrData'. Type '?mederrData' for help.")
	}
	
	if (dat@size > 100000) cat("Preparing the data...")
	data <- dat@data
	err.dim <- dat@numi
	harm.dim <- 2
	hosp.dim <- dat@numj
	
	dat.sort <- data[sort(data$groupj, index.return = TRUE)$ix, ]
	tbl <- array(NA, c(harm.dim, err.dim, hosp.dim))
	for (j in 1:hosp.dim) tbl[, , j] <- t(dat.sort[(1 + err.dim*(j - 1)):(err.dim*j), 1:2])
	tbl[2, , ] <- tbl[2, , ] - tbl[1, , ]
	
	marg.err <- apply(tbl, c(1, 3), sum)
	marg.harm <- apply(tbl, c(2, 3), sum)
	joint.err.harm <- apply(tbl, c(1, 2), sum)
	marg.errharm <- colSums(marg.harm)

	if (stratified) {
		ind.joint.err.harm <- matrix(as.vector(t(marg.err))/rep(marg.errharm, harm.dim), nrow = harm.dim, ncol = hosp.dim, byrow = TRUE) %*% t(marg.harm)
	}
	else {
		marg.err <- colSums(joint.err.harm)
		marg.harm <- rowSums(joint.err.harm)
		total <- sum(joint.err.harm)
		ind.joint.err.harm <- marg.harm%*%t(marg.err)/total
	}

	x <- as.vector(t(joint.err.harm))
	E <- as.vector(t(ind.joint.err.harm))
	if (dat@size > 100000) cat("...done!\n")
	
	n <- length(x)
	error <- 1
	j <- 0
	alpha <- theta0[1]
	beta <- theta0[2]
	theta.eb <- theta0
	
	while((error > toler) & (j < maxiter)) {
		olda <- alpha
		oldb <- beta
		oldll <- negbinom.loglik(c(olda, oldb), x, E)
		oldscore <- negbinom.score(c(alpha, beta), x, E)
		if (j == 0) print(paste("Iteration 1: ||score|| = ", formatC(sqrt(sum(oldscore^2)), digits = 4, format = "f"), " - loglik = ",
							formatC(oldll, digits = 4, format = "f"), sep = ""), quote = FALSE)
		
		alpha <- BBsolve(olda, d.loglik.alpha, N = x, E = E, theta = oldb, quiet = TRUE)$par
		beta <- BBsolve(oldb, d.loglik.beta, N = x, E = E, theta = alpha, quiet = TRUE)$par

		theta.eb <- c(alpha, beta)
		if (any(theta.eb <= 0)) {
			stop("Initial values not valid. Try a different set. Type '?negbinom.em' for help.")
		}

		loglik <- negbinom.loglik(theta.eb, x, E)
		score <- negbinom.score(theta.eb, x, E)

		error <- abs(oldll - loglik)
		j <- j + 1
		if ((j/100 == round(j/100))) print(paste("iter ", j, " / ", maxiter, " - error = ", formatC(error, digits = 5, format = "f"), " > ",
							toler, " - ||score|| = ", formatC(sqrt(sum(score^2)), digits = 4, format = "f"), " - loglik = ",
							formatC(loglik, digits = 4, format = "f"), sep = ""), quote = FALSE)
	}
	print(paste("Iteration ", (j - 1)," (last): ||score|| = ", formatC(sqrt(sum(negbinom.score(theta.eb, x, E)^2)), digits = 4, format = "f"),
					" - loglik = ", formatC(negbinom.loglik(theta.eb, x, E), digits = 4, format = "f"), sep = ""), quote = FALSE)

	if (se) {
		hess.lower <- genD(negbinom.loglik, theta.eb, N = x, E = E)$D[, 3:5]
		hess <- matrix(NA, nrow = 2, ncol = 2)
		for (i in 1:2) {
			hess[i, (1:i)] <- hess.lower[(1 + sum(0:(i - 1))):sum(1:i)]
			hess[(1:i), i] <- hess[i, (1:i)]
		}
		se <- sqrt(diag(solve(-hess)))
		out <- list(theta.hat = theta.eb, final.err = error, final.ll = loglik, final.score = score, num.iter = (j - 1), se = se, N = x, E = E,
				prior = "gamma")
	}
	else {
		out <- list(theta.hat = theta.eb, final.err = error, final.ll = loglik, final.score = score, num.iter = (j - 1), se = NA, N = x, E = E,
				prior = "gamma")
	}

	return(out)
}

negbinom.score <- function(theta, N, E) {
	score <- rep(NA, 2)
	score[1] <- sum(-log(1 + E/theta[2]) + digamma(theta[1] + N) - digamma(theta[1]))
	score[2] <- sum(-1/(E + theta[2])*(N - theta[1]*E/theta[2]))
	return(score)
}
