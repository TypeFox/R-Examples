FDRreg = function(z, covars, nulltype='empirical', type='linear', nmc=10000, nburn=500, nmids=150, densknots=10, regknots=5) {
# False discovery rate regression
# z = vector of z scores
# covars = design matrix of covariates, assumed NOT to have an intercept just as in vanilla lm()
# nulltype = flag for what kind of null hypothesis to assume, theoretical or empirical
	N = length(z)
	fjindex=list()
	if(type=='linear') {
		X = cbind(1, colnormalize(covars))
		fjindex = as.list(2:ncol(X))
	}
	else if(type=='additive') {
		X = matrix(rep(1,N))
		for(j in 1:ncol(covars)) {
			Xj = make.bspline.matrix(as.matrix(covars[,j]), nknots=regknots, method='equi')
			fjindex[[j]] = ncol(X) + 1:ncol(Xj)
			X = cbind(X, Xj)
		}
	}
	else {
		warning('Invalid model type specified: assuming linearity in the covariates.')
	}
	P = ncol(X)

	# Estimate the marginal density and empirical null
	l1 = efron(z, nmids=nmids, df=densknots)
	mu0 = l1$mu0
	sig0 = l1$sig0
	if(nulltype=='empirical') {
		M0 = dnorm(z,mu0,sig0)
		grid.f0z = dnorm(l1$mids,mu0,sig0)
	}
	else if(nulltype=='theoretical') {
		M0 = dnorm(z)
		grid.f0z = dnorm(l1$mids)
	}
	else {
		warning('Invalid null type specified: assuming empirical null.')
		M0 = dnorm(z,mu0,sig0)
		grid.f0z = dnorm(l1$mids,mu0,sig0)
	}
	MTot = l1$fz
	nullcases = which(M0>MTot)

	# Initialize MCMC
	PriorPrec = diag(c(1e-9, rep(0.1, P-1)))
	PriorPrecXMean = rep(0,P)
	
	# Pass to C++ and return result
	out1 = FDRregCPP(z, X, M0=M0, MTot = MTot, PriorPrec, PriorPrecXMean, nmc=nmc, nburn=nburn)
	out2 = getFDR(out1$postprob)
	list(z=z, localfdr=out2$localfdr, FDR=out2$FDR, X=X, grid=l1$mids, breaks=l1$breaks,
         grid.fz=l1$zdens, grid.f0z=grid.f0z, grid.zcounts=l1$zcounts, dnull = M0,
         dmix=MTot, empirical.null=list(mu0=mu0, sig0=sig0), betasave = out1$betasave,
         priorprob = out1$priorprob, postprob = out1$postprob, fjindex=fjindex)
	
	# Beta = rep(0, P)
	# BetaSave = matrix(0, nrow=NMC, ncol=P)
	# PostProbSave = 0
	# PriorProbSave = 0
	
	# # Main MCMC
	# for(t in 1:(NMC+burn)) {
		# # Update indicators
		# Psi = drop(X %*% Beta)
		# W = ilogit(Psi)
		# PostProb = 1-pmin(1, {1-W}*M0/MTot)
		# PostProb[nullcases] = 0
		# Gamma = rbinom(N,1,PostProb)
			
		# # Update latent variables in logit likelihood
		# Om = as.numeric(rpg.devroye(N,rep(1,N),Psi))
		  
		# # Update regression parameters
		# Kap = PostProb - 1/2
		# PrecMat = t(X) %*% {Om * X} + PriorPrec
		# Beta.V = solve(PrecMat)
		# Beta.mu = Beta.V %*% {t(X) %*% Kap + PriorPrecXMean}
		# Beta = t(rmvnorm(1,mean=Beta.mu,sigma=Beta.V))	
		# if(t > burn) {
			# BetaSave[t-burn,] = Beta
			# PostProbSave = PostProbSave + (1/NMC)*PostProb
			# PriorProbSave = PriorProbSave + (1/NMC)*W
		# }
	# }
	# list(mids=l1$mids, zcounts=l1$zcounts, zdens=l1$zdens, fz=l1$fz, mu0=mu0, sig0=sig0,
		# X=X, BetaSave = BetaSave, PriorProb = PriorProbSave, PostProb = PostProbSave, fjindex=fjindex)
}
