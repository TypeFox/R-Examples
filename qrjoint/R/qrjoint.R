## version 2. Playing with the default for kappa

qrjoint <- function(x, y, nsamp = 1e3, thin = 10, cens = rep(0,length(y)), incr = 0.01, par = "prior", nknots = 6, hyper = list(sig = c(.1,.1), lam = c(6,4), kap = c(0.1,0.1,1)), shrink = FALSE, prox.range = c(.2,.95), acpt.target = 0.15, ref.size = 3, blocking = "std5", temp = 1, expo = 2, blocks.mu, blocks.S, fix.nu = FALSE){
	
	x <- as.matrix(x)
    n <- length(y); x <- matrix(x, nrow = n); p <- ncol(x)
	x.names <- dimnames(x)[[2]]; if(is.null(x.names)) x.names <- paste("X", 1:p, sep = "")
	x <- scale(x, chull.center(x))
	
	Ltail <- ceiling(log(n*incr,2))
	if(Ltail > 0){
		tau.tail <- incr / 2^(Ltail:1)
		tau.g <- c(0, tau.tail, seq(incr, 1 - incr, incr), 1 - tau.tail[Ltail:1], 1)
		L <- length(tau.g); mid <- which.min(abs(tau.g - 0.5)); reg.ix <- (1:L)[-c(1 + 1:Ltail, L - 1:Ltail)]
	} else {
		tau.g <- seq(0, 1, incr)
		L <- length(tau.g); mid <- which.min(abs(tau.g - 0.5)); reg.ix <- (1:L)
	}
	
	tau.kb <- seq(0,1,len = nknots) 
	tau.k <- tau.kb
	
	a.sig <- hyper$sig; if(is.null(a.sig)) a.sig <- c(.1, .1)
	a.lam <- hyper$lam; if(is.null(a.lam)) a.lam <- c(6, 4)
	a.kap <- hyper$kap; if(is.null(a.kap)) a.kap <- c(1.5, 1.5, 1); a.kap <- matrix(a.kap, nrow = 3); nkap <- ncol(a.kap); a.kap[3,] <- log(a.kap[3,])
	hyper <- c(a.sig, c(a.kap))
	
	prox.grid <- proxFn(max(prox.range), min(prox.range), 0.5)
	ngrid <- length(prox.grid)
	lamsq.grid <- lamFn(prox.grid)^2
	prior.grid <- -diff(pbeta(c(1, (prox.grid[-1] + prox.grid[-ngrid])/2, 0), a.lam[1], a.lam[2]))
	lp.grid <- log(prior.grid)
	
	d.kg <- abs(outer(tau.k, tau.g, "-"))^expo	
	d.kk <- abs(outer(tau.k, tau.k, "-"))^expo
	gridmats <- matrix(NA, nknots*(L + nknots)+2, ngrid)
	K0 <- 0
    t1 <- Sys.time()
	for(i in 1:ngrid){
		K.grid <- exp(-lamsq.grid[i] * d.kg); K.knot <- exp(-lamsq.grid[i] * d.kk);	diag(K.knot) <- 1 + 1e-10	
		R.knot <- chol(K.knot); A.knot <- solve(K.knot, K.grid)
		gridmats[,i] <- c(c(A.knot), c(R.knot), sum(log(diag(R.knot))), lp.grid[i])		
		K0 <- K0 + prior.grid[i] * K.knot
	}
    t2 <- Sys.time()
    cat("Matrix calculation time per 1e3 iterations =", round(1e3 * as.numeric(t2 - t1), 2), "\n")
	
	niter <- nsamp * thin
	dimpars <- c(n, p, L, mid - 1, nknots, ngrid, ncol(a.kap), niter, thin, nsamp)
	
	if(par[1] == "prior") {
		par <- rep(0, (nknots+1) * (p+1) + 2)
		if(fix.nu) par[(nknots+1) * (p+1) + 2] <- nuFn.inv(fix.nu)
		
		beta.rq <- sapply(tau.g, function(a) return(coef(rq(y ~ x, tau = a))))
		v <- bs(tau.g, df = 5)
		rq.lm <- apply(beta.rq, 1, function(z) return(coef(lm(z ~ v))))
		
		delta <- tau.g[2]
		tau.0 <- tau.g[mid]
		rq.tau0 <- c(c(1, predict(v, tau.0)) %*% rq.lm)
		rq.delta <- c(c(1, predict(v, delta)) %*% rq.lm)
		rq.deltac <- c(c(1, predict(v, 1 - delta)) %*% rq.lm)
		
		par[nknots*(p+1) + 1:(p+1)] <- as.numeric(rq.tau0)
		sigma <- 1
		par[(nknots+1)*(p+1) + 1] <- sigFn.inv(sigma, a.sig)
		
		for(i in 1:(p+1)){
			kapsq <- sum(exp(a.kap[3,]) * (a.kap[2,] / a.kap[1,]))
			lam.ix <- sample(length(lamsq.grid), 1, prob = prior.grid)
			R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),lam.ix], nknots, nknots)
			z <- sqrt(kapsq) * c(crossprod(R, rnorm(nknots)))
			par[(i - 1) * nknots + 1:nknots] <- z - mean(z)
			
		}
		
		beta.hat <- estFn(par, x, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, FALSE)
		qhat <- tcrossprod(cbind(1, x), beta.hat)
		
		infl <- max(max((y - qhat[,mid])/(qhat[,ncol(qhat) - 1] - qhat[,mid])), max((qhat[,mid] - y)/(qhat[,mid] - qhat[,2])))
		oo <- .C("INIT", par = as.double(par), x = as.double(x), y = as.double(y), cens = as.integer(cens),
				 shrink = as.integer(shrink), hyper = as.double(hyper), dim = as.integer(dimpars), gridpars = as.double(gridmats), 
				 tau.g = as.double(tau.g), siglim = as.double(sigFn.inv(c(1.0 * infl * sigma, 10 * infl * sigma), a.sig)))
		
		par <- oo$par		
	} else if (par[1] == "RQ"){
		
		par <- rep(0, (nknots+1) * (p+1) + 2)
		
		beta.rq <- sapply(tau.g, function(a) return(coef(rq(y ~ x, tau = a))))
		v <- bs(tau.g, df = 5)
		rq.lm <- apply(beta.rq, 1, function(z) return(coef(lm(z ~ v))))
		
		delta <- tau.g[2]
		tau.0 <- tau.g[mid]
		rq.tau0 <- c(c(1, predict(v, tau.0)) %*% rq.lm)
		rq.delta <- c(c(1, predict(v, delta)) %*% rq.lm)
		rq.deltac <- c(c(1, predict(v, 1 - delta)) %*% rq.lm)
		
		par[nknots*(p+1) + 1:(p+1)] <- as.numeric(rq.tau0)
		nu <- ifelse(fix.nu, fix.nu, nuFn(0))
		sigma <- min((rq.delta[1] - rq.tau0[1]) / Q0(delta, nu), (rq.deltac[1] - rq.tau0[1]) / Q0(1 - delta, nu))
		par[(nknots+1)*(p+1) + 1]  <- sigFn.inv(sigma, a.sig)
		
		epsilon <- 0.1 * min(diff(sort(tau.k)))
		tau.knot.plus <- pmin(tau.k + epsilon, 1)
		tau.knot.minus <- pmax(tau.k - epsilon, 0)
		beta.rq.plus <- cbind(1, predict(v, tau.knot.plus)) %*% rq.lm
		beta.rq.minus <- cbind(1, predict(v, tau.knot.minus)) %*% rq.lm
		zeta0.plus <- F0((beta.rq.plus[,1] - rq.tau0[1]) / sigma, nu)
		zeta0.minus <- F0((beta.rq.minus[,1] - rq.tau0[1]) / sigma, nu)
		zeta0.dot.knot <- (zeta0.plus - zeta0.minus) / (tau.knot.plus - tau.knot.minus)
		w0.knot <- log(pmax(epsilon, zeta0.dot.knot)) / shrinkFn(p)
		w0.knot <- (w0.knot - mean(w0.knot)) 
		
		w0PP <- ppFn0(w0.knot, gridmats, L, nknots, ngrid)
		w0 <- w0PP$w	
		
		zeta0.dot <- exp(shrinkFn(p) * (w0 - max(w0)))
		zeta0 <- trape(zeta0.dot[-c(1,L)], tau.g[-c(1,L)], L-2)
		zeta0.tot <- zeta0[L-2] 
		zeta0 <- c(0, tau.g[2] + (tau.g[L-1]-tau.g[2])*zeta0 / zeta0.tot, 1)
		zeta0.dot <- (tau.g[L-1]-tau.g[2])*zeta0.dot / zeta0.tot
		zeta0.dot[c(1,L)] <- 0
		
		par[1:nknots] <- w0.knot		
		beta0.dot <- sigma * q0(zeta0, nu) * zeta0.dot
		
		tilt.knot <- tau.g[tilt.ix <- sapply(tau.k, function(a) which.min(abs(a - zeta0)))]
		
		tau.knot.plus <- pmin(tilt.knot + epsilon, 1)
		tau.knot.minus <- pmax(tilt.knot - epsilon, 0)
		beta.rq.plus <- cbind(1, predict(v, tau.knot.plus)) %*% rq.lm
		beta.rq.minus <- cbind(1, predict(v, tau.knot.minus)) %*% rq.lm
		
		beta.dot.knot <- (beta.rq.plus[,-1,drop=FALSE] - beta.rq.minus[,-1,drop=FALSE]) /  (tau.knot.plus - tau.knot.minus)
		
		par[nknots + 1:(nknots*p)] <- c(beta.dot.knot)
		beta.hat <- estFn(par, x, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, FALSE)
		qhat <- tcrossprod(cbind(1, x), beta.hat)
		infl <- max(max((y - qhat[,mid])/(qhat[,ncol(qhat) - 1] - qhat[,mid])),
					max((qhat[,mid] - y)/(qhat[,mid] - qhat[,2])))
		
		oo <- .C("INIT", par = as.double(par), x = as.double(x), y = as.double(y), cens = as.integer(cens),
				 shrink = as.integer(shrink), hyper = as.double(hyper), dim = as.integer(dimpars), gridpars = as.double(gridmats), 
				 tau.g = as.double(tau.g), siglim = as.double(sigFn.inv(c(1.0 * infl * sigma, 10.0 * infl * sigma), a.sig)))
		
		par <- oo$par		
	}	
	
	npar <- (nknots+1)*(p+1) + 2
	if(blocking == "single"){
		blocks <- list(rep(TRUE, npar))
	} else if(blocking == "single2"){
		blocks <- list(rep(TRUE, npar), rep(FALSE, npar))
		blocks[[2]][nknots*(p+1) + 1:(p+3)] <- TRUE
	} else if(blocking == "single3"){
		blocks <- list(rep(TRUE, npar), rep(FALSE, npar), rep(FALSE, npar))
		blocks[[2]][nknots*(p+1) + 1:(p+1)] <- TRUE
		blocks[[3]][(nknots+1)*(p+1) + 1:2] <- TRUE
	} else if(blocking == "std0"){
		blocks <- replicate(p + 1, rep(FALSE, npar), simplify = FALSE)
		for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
	} else if(blocking == "std1"){
		blocks <- replicate(p + 2, rep(FALSE, npar), simplify = FALSE)
		for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
		blocks[[p + 2]][nknots*(p+1) + 1:(p+3)] <- TRUE
	} else if(blocking == "std2"){
		blocks <- replicate(p + 3, rep(FALSE, npar), simplify = FALSE)
		for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1, (nknots+1)*(p+1) + 1:2)] <- TRUE
		blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE
		blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE
	} else if(blocking == "std3"){
		blocks <- replicate(p + 3, rep(FALSE, npar), simplify = FALSE)
		for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots)] <- TRUE
		blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE
		blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE
	} else if(blocking == "std4"){
		blocks <- replicate(p + 3, rep(FALSE, npar), simplify = FALSE)
		for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1)] <- TRUE
		blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE
		blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE
	} else if(blocking == "std5"){
		blocks <- replicate(p + 4, rep(FALSE, npar), simplify = FALSE)
		for(i in 0:p) blocks[[i + 1]][c(i * nknots + 1:nknots, nknots*(p+1) + i + 1)] <- TRUE
		blocks[[p + 2]][nknots*(p+1) + 1:(p+1)] <- TRUE
		blocks[[p + 3]][(nknots+1)*(p+1) + 1:2] <- TRUE
		blocks[[p+4]][1:npar] <- TRUE
	} else {
		blocks <- replicate(npar, rep(FALSE, npar), simplify = FALSE)
		for(i in 1:npar) blocks[[i]][i] <- TRUE
	}
	
	
	nblocks <- length(blocks)
	if(fix.nu) for(j in 1:nblocks) blocks[[j]][(nknots+1) * (p+1) + 2] <- FALSE
	
	blocks.ix <- c(unlist(lapply(blocks, which))) - 1
	blocks.size <- sapply(blocks, sum)
	if(missing(blocks.mu)) blocks.mu <- rep(0, sum(blocks.size))
	if(missing(blocks.S)){
		blocks.S <- lapply(blocks.size, function(q) diag(1, q))
		if(substr(blocking, 1, 4) == "std"){
			for(i in 1:(p+1)) blocks.S[[i]][1:nknots, 1:nknots] <- K0
			if(as.numeric(substr(blocking, 5,6)) > 1){
				blocks.S[[p + 2]] <- summary(rq(y ~ x, tau = 0.5), se = "boot", cov = TRUE)$cov
				blocks.S[[p + 3]] <- matrix(c(1, 0, 0, .1), 2, 2)
			}
			if(as.numeric(substr(blocking, 5,6)) == 5){
				slist <- list(); length(slist) <- p + 3
				for(i in 1:(p+1)) slist[[i]] <- K0
				slist[[p+2]] <- summary(rq(y ~ x, tau = 0.5), se = "boot", cov = TRUE)$cov
				slist[[p+3]] <- matrix(c(1, 0, 0, .1), 2, 2)
				blocks.S[[p + 4]] <- as.matrix(bdiag(slist))
			}
		}
		blocks.S <- unlist(blocks.S)
	}
		
	imcmc.par <- c(nblocks, ref.size, TRUE, max(10, niter/1e4), rep(0, nblocks))
	dmcmc.par <- c(temp, 0.999, rep(acpt.target, nblocks), 2.38 / sqrt(blocks.size))
	
	
	tm.c <- system.time(oo <- .C("BJQR", par = as.double(par), x = as.double(x), y = as.double(y), cens = as.integer(cens),
								 shrink = as.integer(shrink), hyper = as.double(hyper), dim = as.integer(dimpars), gridmats = as.double(gridmats), 
								 tau.g = as.double(tau.g), muV = as.double(blocks.mu), SV = as.double(blocks.S), blocks = as.integer(blocks.ix), 
								 blocks.size = as.integer(blocks.size), dmcmcpar = as.double(dmcmc.par), 
								 imcmcpar = as.integer(imcmc.par), parsamp = double(nsamp * length(par)), 
								 acptsamp = double(nsamp * nblocks), lpsamp = double(nsamp)))  
	cat("elapsed time:", round(tm.c[3]), "seconds\n")
	
	oo$x <- x; oo$y <- y; oo$xnames <- x.names; oo$gridmats <- gridmats; oo$prox <- prox.grid; oo$reg.ix <- reg.ix; oo$runtime <- tm.c[3]
	class(oo) <- "qrjoint"
	return(oo)
}

update.qrjoint <- function(object, nadd, append = TRUE, ...){
	niter <- object$dim[8]; thin <- object$dim[9]; nsamp <- object$dim[10]
	if(missing(nadd)) nadd <- nsamp
	par <- object$par; npar <- length(par)
	dimpars <- object$dim
	dimpars[8] <- nadd * thin
	dimpars[10] <- nadd
	nblocks <- object$imcmcpar[1]
	object$imcmcpar[4] <- max(10, nadd * thin/1e4)
	
	tm.c <- system.time(oo <- .C("BJQR", par = as.double(par), x = as.double(object$x), y = as.double(object$y), cens = as.integer(object$cens),
								 shrink = as.integer(object$shrink), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), 
								 tau.g = as.double(object$tau.g), muV = as.double(object$muV), SV = as.double(object$SV), blocks = as.integer(object$blocks), 
								 blocks.size = as.integer(object$blocks.size), dmcmcpar = as.double(object$dmcmcpar), 
								 imcmcpar = as.integer(object$imcmcpar), parsamp = double(nadd * npar), 
								 acptsamp = double(nadd * nblocks), lpsamp = double(nadd)))  
	cat("elapsed time:", round(tm.c[3]), "seconds\n")
	
	oo$x <- object$x; oo$y <- object$y; oo$xnames <- object$xnames; oo$gridmats <- object$gridmats; oo$prox <- object$prox; oo$reg.ix <- object$reg.ix; oo$runtime <- object$runtime+tm.c[3]
	if(append){
		oo$dim[8] <- niter + nadd * thin
		oo$dim[10] <- nsamp + nadd
		oo$parsamp <- c(object$parsamp, oo$parsamp)
		oo$acptsamp <- c(object$acptsamp, oo$acptsamp)
		oo$lpsamp <- c(object$lpsamp, oo$lpsamp)
	}
	class(oo) <- "qrjoint"
	return(oo)		
}

coef.qrjoint <- function(object, burn.perc = 0.5, nmc = 200, plot = FALSE, show.intercept = TRUE, reduce = TRUE, ...){
	niter <- object$dim[8]	
	nsamp <- object$dim[10]
	pars <- matrix(object$parsamp, ncol = nsamp)
	ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 1)[-1]))
	
	n <- object$dim[1]; p <- object$dim[2]; L <- object$dim[3]; mid <- object$dim[4] + 1; nknots <- object$dim[5]; ngrid <- object$dim[6]
	a.sig <- object$hyper[1:2]; a.kap <- matrix(object$hyper[-c(1:2)], nrow = 3)
	tau.g <- object$tau.g; reg.ix <- object$reg.ix
	x.ce <- outer(rep(1, L), attr(object$x, "scaled:center")); x.sc <- outer(rep(1,L), attr(object$x, "scaled:scale"))

	beta.samp <- apply(pars[,ss], 2, function(p1) c(estFn(p1, object$x, object$y, object$gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, reduce, x.ce, x.sc)))

	if(reduce) tau.g <- tau.g[reg.ix]
	L <- length(tau.g)
	if(plot){
		nr <- ceiling(sqrt(p+show.intercept)); nc <- ceiling((p+show.intercept)/nr)
		par(mfrow = c(nr, nc))
	}
	reach <- 0
	beta.hat <- list()
	plot.titles <- c("Intercept", object$xnames)
	j <- 1
	b <- beta.samp[reach + 1:L,]
	beta.hat[[j]] <- getBands(b, plot = (plot & show.intercept), add = FALSE, x = tau.g, xlab = "tau", ylab = "Coefficient", bty = "n", ...)
	if(plot & show.intercept) title(main = plot.titles[j])
	reach <- reach + L
	for(j in 2:(p+1)){
		b <- beta.samp[reach + 1:L,]
		beta.hat[[j]] <- getBands(b, plot = plot, add = FALSE, x = tau.g, xlab = "tau", ylab = "Coefficient", bty = "n", ...)
		if(plot) {
			title(main = plot.titles[j])
			abline(h = 0, lty = 2, col = 4)
		}
		reach <- reach + L
	}	
	names(beta.hat) <- plot.titles
	invisible(list(beta.samp = beta.samp, beta.est = beta.hat))
}

summary.qrjoint <- function(object, ntrace = 1000, plot.dev = TRUE, more.details = FALSE, ...){
	thin <- object$dim[9]	
	nsamp <- object$dim[10]
	pars <- matrix(object$parsamp, ncol = nsamp)
	ss <- unique(pmax(1, round(nsamp * (1:ntrace/ntrace))))
	dimpars <- object$dim
	dimpars[8] <- length(ss)
	
	n <- object$dim[1]; p <- object$dim[2]; ngrid <- object$dim[6]
	sm <- .C("DEV", pars = as.double(pars[,ss]), x = as.double(object$x), y = as.double(object$y), cens = as.integer(object$cens),
			 shrink = as.integer(object$shrink), hyper = as.double(object$hyper), dim = as.integer(dimpars), gridmats = as.double(object$gridmats), tau.g = as.double(object$tau.g),
			 devsamp = double(length(ss)), llsamp = double(length(ss)*n), pgsamp = double(length(ss)*ngrid*(p+1)))
	deviance <- sm$devsamp
	ll <- matrix(sm$llsamp, ncol = length(ss))
	fit.waic <- waic(ll)
	pg <- matrix(sm$pgsamp, ncol = length(ss))
	prox.samp <- matrix(NA, p+1, length(ss))
	for(i in 1:(p+1)){
		prox.samp[i,] <- object$prox[apply(pg[(i-1)*ngrid + 1:ngrid,], 2, function(pr) sample(length(pr), 1, prob = pr))]
	}	
	
	if(more.details) par(mfrow = c(2,2), mar = c(5,4,3,2)+.1)
	if(plot.dev){
		plot(thin * ss, deviance, ty = "l", xlab = "Markov chain iteration", ylab = "Deviance", bty = "n", main = "Fit trace plot", ...)
		grid(col = "gray")
	}
	
	if(more.details){
		ngrid <- length(object$prox)
		prior.grid <- exp(object$gridmats[nrow(object$gridmats),])
		lam.priorsamp <- lamFn(sample(object$prox, ntrace, replace = TRUE, prob = prior.grid))
		lam.prior.q <- quantile(lam.priorsamp, pr = c(.025, .5, .975))
		lam.samp <- lamFn(prox.samp)
		a <- min(lamFn(object$prox))
		b <- diff(range(lamFn(object$prox))) * 1.2
		plot(thin * ss, lam.samp[1,], ty = "n", ylim = a + c(0, b * (p + 1)), bty = "n", ann = FALSE, axes = FALSE)
		axis(1)
		for(i in 1:(p+1)){
			abline(h = b * (i-1) + lamFn(object$prox), col = "gray")
			abline(h = b * (i - 1) + lam.prior.q, col = "red", lty = c(2,1,2))
			lines(thin * ss, b * (i-1) + lam.samp[i,], lwd = 1, col = 4)
			if(i %% 2) axis(2, at = b * (i-1) + lamFn(object$prox[c(1,ngrid)]), labels = round(object$prox[c(1,ngrid)],2), las = 1, cex.axis = 0.6) 
			mtext(substitute(beta[index], list(index = i - 1)), side = 4, line = 0.5, at = a + b * (i - 1) + 0.4*b, las = 1)			
		}
		title(xlab = "Markov chain iteration", ylab = "Proxmity posterior", main = "Mixing over GP scaling")					
		
		theta <- as.mcmc(t(matrix(object$parsamp, ncol = nsamp)[,unique(ceiling(seq(0.5*nsamp,nsamp,1)))]))
		gg <- geweke.diag(theta, .1, .5)
		zvals <- gg$z

		pp <- 2 * (1 - pnorm(abs(zvals)))
		plot(sort(pp), ylab = "Geweke p-values", xlab = "Parameter index (reordered)", main = "Convergence diagnosis", ty = "h", col = 4, ylim = c(0, 0.3), lwd = 2)
		abline(h = 0.05, col = 2, lty = 2)
		abline(a = 0, b = 0.1 / length(pp), col = 2, lty = 2)
		mtext(c("BH-10%", "5%"), side = 4, at = c(0.1, 0.05), line = 0.1, las = 0, cex = 0.6)
		
		npar <- length(object$par)
		image(1:npar, 1:npar, cor(theta), xlab = "Parameter index", ylab = "Parameter index", main = "Parameter correlation")
		
	}
	invisible(list(deviance = deviance, pg = pg, prox = prox.samp, ll = ll, waic = fit.waic))
}

estFn <- function(par, x, y, gridmats, L, mid, nknots, ngrid, a.kap, a.sig, tau.g, reg.ix, reduce = TRUE, x.ce = 0, x.sc = 1){
	
	n <- length(y); p <- ncol(x)
	wKnot <- matrix(par[1:(nknots*(p+1))], nrow = nknots)
	w0PP <- ppFn0(wKnot[,1], gridmats, L, nknots, ngrid)
	w0 <- w0PP$w
	wPP <- apply(wKnot[,-1,drop=FALSE], 2, ppFn, gridmats = gridmats, L = L, nknots = nknots, ngrid = ngrid, a.kap = a.kap)
	vMat <- matrix(sapply(wPP, extract, vn = "w"), ncol = p)
	
	zeta0.dot <- exp(shrinkFn(p) * (w0 - max(w0)))
	zeta0 <- trape(zeta0.dot[-c(1,L)], tau.g[-c(1,L)], L-2)
	zeta0.tot <- zeta0[L-2] 
	zeta0 <- c(0, tau.g[2] + (tau.g[L-1]-tau.g[2])*zeta0 / zeta0.tot, 1)
	zeta0.dot <- (tau.g[L-1]-tau.g[2])*zeta0.dot / zeta0.tot
	zeta0.dot[c(1,L)] <- 0
	
	reach <- nknots*(p+1)
	gam0 <- par[reach + 1]; reach <- reach + 1
	gam <- par[reach + 1:p]; reach <- reach + p
	sigma <- sigFn(par[reach + 1], a.sig); reach <- reach + 1
	nu <- nuFn(par[reach + 1]);
	
	b0dot <- sigma * q0(zeta0, nu) * zeta0.dot
	beta0.hat <- rep(NA, L)
	beta0.hat[mid:L] <- gam0 + trape(b0dot[mid:L], tau.g[mid:L], L - mid + 1)
	beta0.hat[mid:1] <- gam0 + trape(b0dot[mid:1], tau.g[mid:1], mid)
	
	vNorm <- sqrt(rowSums(vMat^2))
	a <- tcrossprod(vMat, x)
	aX <- apply(-a, 1, max)/vNorm
	aX[is.nan(aX)] <- Inf
	aTilde <- vMat / (aX * sqrt(1 + vNorm^2))
	ab0 <- b0dot * aTilde
	
	beta.hat <- kronecker(rep(1,L), t(gam))
	beta.hat[mid:L,] <- beta.hat[mid:L,] + apply(ab0[mid:L,,drop=FALSE], 2, trape, h = tau.g[mid:L], len = L - mid + 1)
	beta.hat[mid:1,] <- beta.hat[mid:1,] + apply(ab0[mid:1,,drop=FALSE], 2, trape, h = tau.g[mid:1], len = mid)
	beta.hat <- beta.hat / x.sc
	beta0.hat <- beta0.hat - rowSums(beta.hat * x.ce)
	betas <- cbind(beta0.hat, beta.hat)
	if(reduce) betas <- betas[reg.ix,,drop = FALSE]
	return(betas)
}


chull.center <- function (x, maxEPts = ncol(x) + 1, plot = FALSE){
    sx <- as.matrix(apply(x, 2, function(s) punif(s, min(s), max(s))))
    dd <- rowSums(scale(sx)^2)
    ix.dd <- order(dd, decreasing = TRUE)
    sx <- sx[ix.dd, , drop = FALSE]
    x.chol <- inchol(sx, maxiter = maxEPts)
    ix.epts <- ix.dd[pivots(x.chol)]
    x.choose <- x[ix.epts, , drop = FALSE]
    xCent <- as.numeric(colMeans(x.choose))
    attr(xCent, "EPts") <- ix.epts
    if (plot) {
        n <- nrow(x)
        p <- ncol(x)
        xjit <- x + matrix(rnorm(n * p), n, p) %*% diag(0.05 * apply(x, 2, sd), p)
        xmean <- colMeans(x)
        x.ept <- x[ix.epts, ]
        M <- choose(p, 2)
        xnames <- dimnames(x)[[2]]
        if (is.null(xnames)) xnames <- paste("x", 1:p, sep = "")
        par(mfrow = c(ceiling(M/ceiling(sqrt(M))), ceiling(sqrt(M))), mar = c(2, 3, 0, 0) + 0.1)
        xmax <- apply(x, 2, max)
        for (i in 1:(p - 1)) {
            for (j in (i + 1):p) {
                plot(x[, i], x[, j], col = "gray", cex = 1, ann = FALSE, ty = "n", axes = FALSE, bty = "l")
                title(xlab = xnames[i], line = 0.3)
                title(ylab = xnames[j], line = 0.3)
                ept <- chull(x[, i], x[, j])
                polygon(x[ept, i], x[ept, j], col = gray(0.9), border = "white")
                points(xjit[, i], xjit[, j], pch = ".", col = gray(0.6))
                points(xmean[i], xmean[j], col = gray(0), pch = 17, cex = 1)
                points(xCent[i], xCent[j], col = gray(0), pch = 1, cex = 2)
                points(x.ept[, i], x.ept[, j], col = gray(.3), pch = 10, cex = 1.5)
            }
        }
    }
    return(xCent)
}


waic <- function(logliks, print = TRUE){
	lppd <- sum(apply(logliks, 1, logmean))
	p.waic.1 <- 2 * lppd - 2 * sum(apply(logliks, 1, mean))
	p.waic.2 <- sum(apply(logliks, 1, var))
	waic.1 <- -2 * lppd + 2 * p.waic.1
	waic.2 <- -2 * lppd + 2 * p.waic.2
	if(print) cat("WAIC.1 =", round(waic.1, 2), ", WAIC.2 =", round(waic.2, 2), "\n")
	invisible(c(WAIC1 = waic.1, WAIC2 = waic.2))
}



ppFn0 <- function(w.knot, gridmats, L, nknots, ngrid){
	w.grid <- matrix(NA, L, ngrid)
	lpost.grid <- rep(NA, ngrid)
	for(i in 1:ngrid){
		A <- matrix(gridmats[1:(L*nknots),i], nrow = nknots)
		R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),i], nrow = nknots)
		r <- sum.sq(backsolve(R, w.knot, transpose = TRUE))
		w.grid[,i] <- colSums(A * w.knot)
		lpost.grid[i] <- -(0.5*nknots+1.5)*log1p(0.5*r/1.5) - gridmats[nknots*(L+nknots)+1,i] + gridmats[nknots*(L+nknots)+2,i]		
	}
	lpost.sum <- logsum(lpost.grid)
	post.grid <- exp(lpost.grid - lpost.sum)
	w <- c(w.grid %*% post.grid)
	return(list(w = w, lpost.sum = lpost.sum))
}


ppFn <- function(w.knot, gridmats, L, nknots, ngrid, a.kap){
	w.grid <- matrix(NA, L, ngrid)
	lpost.grid <- rep(NA, ngrid)
	for(i in 1:ngrid){
		A <- matrix(gridmats[1:(L*nknots),i], nrow = nknots)
		R <- matrix(gridmats[L*nknots + 1:(nknots*nknots),i], nrow = nknots)
		r <- sum.sq(backsolve(R, w.knot, transpose = TRUE))
		w.grid[,i] <- colSums(A * w.knot)
		lpost.grid[i] <- (logsum(-(nknots/2+a.kap[1,])*log1p(0.5*r/ a.kap[2,]) + a.kap[3,] + lgamma(a.kap[1,]+nknots/2)-lgamma(a.kap[1,])-.5*nknots*log(a.kap[2,]))
						  - gridmats[nknots*(L+nknots)+1,i] + gridmats[nknots*(L+nknots)+2,i])		
	}
	lpost.sum <- logsum(lpost.grid)
	post.grid <- exp(lpost.grid - lpost.sum)
	w <- c(w.grid %*% post.grid)
	return(list(w = w, lpost.sum = lpost.sum))
}


lamFn <- function(prox) return(sqrt(-100*log(prox)))
nuFn <- function(z) return(0.5 + 5.5*exp(z/2)) 
nuFn.inv <- function(nu) return(2*log((nu - 0.5)/5.5))
sigFn <- function(z, a.sig) return(exp(z/2)) 
sigFn.inv <- function(s, a.sig) return(2 * log(s))
unitFn <- function(u) return(pmin(1 - 1e-10, pmax(1e-10, u)))

q0 <- function(u, nu = Inf) return(1 / (dt(qt(unitFn(u), df = nu), df = nu) * qt(.9, df = nu)))
Q0 <- function(u, nu = Inf) return(qt(unitFn(u), df = nu) / qt(.9, df = nu))
F0 <- function(x, nu = Inf) return(pt(x*qt(.9, df = nu), df = nu))

sum.sq <- function(x) return(sum(x^2))
extract <- function(lo, vn) return(lo[[vn]])
logmean <- function(lx) return(max(lx) + log(mean(exp(lx - max(lx)))))
logsum <- function(lx) return(logmean(lx) + log(length(lx)))
shrinkFn <- function(x) return(1) ##(1/(1 + log(x)))
trape <- function(x, h, len = length(x)) return(c(0, cumsum(.5 * (x[-1] + x[-len]) * (h[-1] - h[-len]))))

getBands <- function(b, col = 2, lwd = 1, plot = TRUE, add = FALSE, x = seq(0,1,len=nrow(b)), remove.edges = TRUE, ...){
	colRGB <- col2rgb(col)/255
	colTrans <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha = 0.2)
	b.med <- apply(b, 1, quantile, pr = .5)
	b.lo <- apply(b, 1, quantile, pr = .025)
	b.hi <- apply(b, 1, quantile, pr = 1 - .025)
	L <- nrow(b)
	ss <- 1:L; ss.rev <- L:1
	if(remove.edges){
		ss <- 2:(L-1); ss.rev <- (L-1):2
	}
	if(plot){
		if(!add) plot(x[ss], b.med[ss], ty = "n", ylim = range(c(b.lo[ss], b.hi[ss])), ...)  
		polygon(x[c(ss, ss.rev)], c(b.lo[ss], b.hi[ss.rev]), col = colTrans, border = colTrans)
		lines(x[ss], b.med[ss], col = col, lwd = lwd)  
	}
	invisible(cbind(b.lo, b.med, b.hi))
}


klGP <- function(lam1, lam2, nknots = 11){
	tau <- seq(0, 1, len = nknots)
	dd <- outer(tau, tau, "-")^2
	K1 <- exp(-lam1^2 * dd); diag(K1) <- 1 + 1e-10; R1 <- chol(K1); log.detR1 <- sum(log(diag(R1)))
	K2 <- exp(-lam2^2 * dd); diag(K2) <- 1 + 1e-10; R2 <- chol(K2); log.detR2 <- sum(log(diag(R2)))
	return(log.detR2-log.detR1 - 0.5 * (nknots - sum(diag(solve(K2, K1)))))
}



proxFn <- function(prox.Max, prox.Min, kl.step = 1){
	prox.grid <- prox.Max
	j <- 1
	while(prox.grid[j] > prox.Min){
		prox1 <- prox.grid[j]
		prox2 <- prox.Min
		kk <- klGP(lamFn(prox1), lamFn(prox2))
		while(kk > kl.step){
			prox2 <- (prox1 + prox2)/2
			kk <- klGP(lamFn(prox1), lamFn(prox2))
		}
		j <- j + 1
		prox.grid <- c(prox.grid, prox2)
	}
	return(prox.grid)
}

