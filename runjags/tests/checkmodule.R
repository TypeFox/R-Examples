library("runjags")
runjags.options(nodata.warning=FALSE)

# Checks that the runjags module distributions are correct.
	
# Loading the dynlib requires rjags for Windows:
if(.Platform$OS.type != 'windows' || requireNamespace('rjags')){

	loaded <- runjags:::dynloadmodule()	
	if(!loaded)
		stop("The internal JAGS dynlib could not be loaded - if you installed this package from CRAN, please file a bug report to the package author")

	cat('Running module tests\n')

	# Required for nchain etc:
	library("coda")

	load("moduletargets.Rsave")

	checksok <- rep(TRUE, N)
	for(i in 1:N){
		success <- try({
		obs <- runjags:::userunjagsmodule(tests[[i]]$distribution, tests[[i]]$funtype, tests[[i]]$parameters, tests[[i]]$x, tests[[i]]$uselog, tests[[i]]$lower)
		expect <- results[[i]]

		# Allow a bit of a larger tolerance than usual as we have saved the results from a different machine:
		problem <- abs(expect-obs) > max(10^-5, .Machine$double.eps^0.5)

		# Or obs is NA:
		problem[is.na(problem)] <- TRUE
		if(any(problem)){
			cat("Error with check number ", i, " (which expected observed):  ", paste(paste(which(problem), " ", expect[problem], " ", obs[problem], ";  ", sep=""), collapse=""), "\n", sep="")
			checksok[i] <- FALSE
		}
		})
		if(inherits(success, 'try-error')){
			checksok[i] <- FALSE
			cat("Uncaught crash error with check number", i, "\n")
		}
	}

	if(!all(checksok)) stop(paste("The runjags module checks failed for test number(s) ", paste(which(!checksok),collapse=","), sep=""))
	
	# From:  Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models. Bayesian Analysis, (3), 515–533.
	# ... the half cauchy should be the same as a ratio of two gammas:
	scale <- seq(1,100,length.out=10)
	quantiles <- seq(0.05,0.95,length.out=10)
	checksok <- rep(TRUE, 10)
	set.seed(1)

	for(i in 1:10){
		success <- try({

		xi <- rnorm(10^6, 0, sd=scale[i])
		tau <- rgamma(10^6, shape=0.5, rate=0.5)  	# Chi-squared with 1 df - i.e.= rchisq(10^6, 1)
		gelman <- abs(xi)/sqrt(tau)
		expect <- quantile(gelman,quantiles)
		obs <- runjags:::userunjagsmodule('halfcauchy', 'q', scale[i], quantiles, FALSE, TRUE)
		# Allow a fairly large 5% tolerance as this is a numerical approximation:
		problem <- (abs(expect-obs)/expect) > 0.05

		# Or obs is NA:
		problem[is.na(problem)] <- TRUE
		if(any(problem)){
			cat("Error with Half Cauchy check number ", i, " - scale=", scale[i], " - (which expected observed):  ", paste(paste(which(problem), " ", expect[problem], " ", obs[problem], ";  ", sep=""), collapse=""), "\n", sep="")
			checksok[i] <- FALSE
		}
		})
		if(inherits(success, 'try-error')){
			checksok[i] <- FALSE
			cat("Uncaught crash error with Half Cauchy check number", i, "\n")
		}		
	}
	if(!all(checksok)) stop(paste("The runjags module Half Cauchy checks failed for test number(s) ", paste(which(!checksok),collapse=","), sep=""))
	
	cat("The internal module tests were passed\n")
}else{
	cat("The internal module tests were skipped (not available on Windows)\n")
}


# Require the rjags library to run the rest of the checks:
dotests <- TRUE
if(!require("rjags")){
	cat("The module checks were not performed as rjags is not installed\n")
	dotests <- FALSE
}

if(dotests){
	
	# Try to load the module:
	loaded <- load.runjagsmodule(fail=FALSE)
	if(!loaded){		
		warning("The module checks were not performed as the internal JAGS module could not be loaded")
		dotests <- FALSE
	}

	# Now check the JAGS implementations just to make sure the functions are found OK:

	m <- "model{
	
		r0 ~ dpar(1,1)
		f0 ~ dpar(1,1)
		d0 <- dpar(0.5,1,1)
		p0 <- ppar(0.5,1,1)
		q0 <- qpar(0.5,1,1)
	
		r1 ~ dpar1(1,1)
		f1 ~ dpar1(1,1)
		d1 <- dpar1(0.5,1,1)
		p1 <- ppar1(0.5,1,1)
		q1 <- qpar1(0.5,1,1)
	
		r2 ~ dpar2(1,1,0)
		f2 ~ dpar2(1,1,0)
		d2 <- dpar2(0.5,1,1,0)
		p2 <- ppar2(0.5,1,1,0)
		q2 <- qpar2(0.5,1,1,0)

		r3 ~ dpar3(1,0,1)
		f3 ~ dpar3(1,0,1)
		d3 <- dpar3(0.5,1,0,1)
		p3 <- ppar3(0.5,1,0,1)
		q3 <- qpar3(0.5,1,0,1)

		r4 ~ dpar4(1,1,0,1)
		f4 ~ dpar4(1,1,0,1)
		d4 <- dpar4(0.5,1,1,0,1)
		p4 <- ppar4(0.5,1,1,0,1)
		q4 <- qpar4(0.5,1,1,0,1)

		rl ~ dlomax(1,1)
		fl ~ dlomax(1,1)
		dl <- dlomax(0.5,1,1)
		pl <- plomax(0.5,1,1)
		ql <- qlomax(0.5,1,1)
	
		rm ~ dmouch(1)
		fm ~ dmouch(1)
		dm <- dmouch(0.5,1)
		pm <- pmouch(0.5,1)
		qm <- qmouch(0.5,1)
	
		rg ~ dgenpar(1,1,1)
		fg ~ dgenpar(1,1,1)
		dg <- dgenpar(0.5,1,1,1)
		pg <- pgenpar(0.5,1,1,1)
		qg <- qgenpar(0.5,1,1,1)

		rh ~ dhalfcauchy(25)
		fh ~ dhalfcauchy(25)
		dh <- dhalfcauchy(25, 25) 
		ph <- phalfcauchy(25, 25) 
		qh <- qhalfcauchy(0.5, 25) 
	
		#monitor# r0,d0,p0,q0,r1,d1,p1,q1,r2,d2,p2,q2,r3,d3,p3,q3,r4,d4,p4,q4,rl,dl,pl,ql,rg,dg,pg,qg,rm,dm,pm,qm,rh,fh,dh,ph,qh
		#data# f0, f1, f2, f3, f4, fl, fg, fm, fh
		#inits# rg
	}"
	rg <- list(1,2)
	f0=f1=f2=f3=f4=fl=fg=fm=fh <- 2

	r <- run.jags(m, n.chains=2, burnin=100, sample=100, method='rjags')
	stopifnot(nchain(as.mcmc.list(r))==2)


	# Now check that my p and q functions are symmetric:

	m <- "model{
	
		qsn <- qnorm(sp, cont1, pos1)
		psn <- pnorm(qsn, cont1, pos1)

		qsg <- qgamma(sp, pos1, pos2)
		psg <- pgamma(qsg, pos1, pos2)
	
		qsl <- qlnorm(sp, cont1, pos2)
		psl <- plnorm(qsl, cont1, pos2)

		qsp <- qpar(sp, pos1, pos2)
		psp <- ppar(qsp, pos1, pos2)
	
		qsp1 <- qpar1(sp, pos1, pos2)
		psp1 <- ppar1(qsp1, pos1, pos2)
	
		qsp2 <- qpar2(sp, pos1, pos2, cont1)
		psp2 <- ppar2(qsp2, pos1, pos2, cont1)
	
		qsp3 <- qpar3(sp, pos1, cont1, pos3)
		psp3 <- ppar3(qsp3, pos1, cont1, pos3)
	
		qsp4 <- qpar4(sp, pos1, pos2, cont1, pos3)
		psp4 <- ppar4(qsp4, pos1, pos2, cont1, pos3)
	
		qsgp <- qgenpar(sp, pos1, cont2, cont1)
		psgp <- pgenpar(qsgp, pos1, cont2, cont1)
	
		qslm <- qlomax(sp, pos1, pos2)
		pslm <- plomax(qslm, pos1, pos2)
	
		qsm <- qmouch(sp, pos1)
		psm <- pmouch(qsm, pos1)
	
		qshc <- qhalfcauchy(sp, pos1)
		pshc <- phalfcauchy(qshc, pos1)
	
		dummy ~ dlomax(1,1)
	
		#monitor# psn, psg, psl, psp1, psp2, psp3, psp4, psgp, pslm, psm, qshc, pshc, 

		#  There appears to be a bug in the code for qpar:	
		#psp,
		#data# sp, cont1, cont2, pos1, pos2, pos3
		#inits# dummy
	}"

	dummy <- 2
	sp <- 1:9/10
	set.seed(1)
	cont1 <- runif(1,-10,10)
	cont2 <- runif(1,-10,10)
	pos1 <- rgamma(1,1,1)
	pos2 <- rgamma(1,1,1)
	pos3 <- rgamma(1,1,1)

	r <- run.jags(m, burnin=100, sample=10, n.chains=1, method='rjags', summarise=FALSE)

	sp <- rep(sp,11)
	ps <- combine.mcmc(r,vars='ps')[1,]
	stopifnot(length(sp)==length(ps))
	if(!isTRUE(all.equal(sp,ps))){
		# Reduce stringency of the test - roughly 10^7 different on linux:
		problem <- abs(sp-ps) > max(10^-4, .Machine$double.eps^0.5)
		if(any(problem)) stop(paste("Error with ps/sp check (which expected observed):  ", paste(paste(names(ps)[which(problem)], " ", sp[problem], " ", ps[problem], ";  ", sep=""), collapse=""), sep=""))
	}
	
	cat("The runjags module tests were passed\n")
	
}else{
	cat("The runjags module tests were skipped (rjags not installed)\n")
}

cat("All module checks passed\n")

