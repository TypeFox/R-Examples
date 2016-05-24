spLMPredictJoint <- function(sp.obj, pred.coords, pred.covars, start = 1, 
	end = nrow(sp.obj$p.theta.samples), thin = 1, verbose = TRUE, n.report = 100, 
	noisy = FALSE, method = "eigen")
{
	nsim <- end - start + 1
	n <- nrow(sp.obj$coords)
	np <- nrow(pred.coords)
	opos <- 1:n
	ppos <- (n+1):(n+np)

	rec <- spBayes::spRecover(sp.obj, get.beta = TRUE, get.w = FALSE, start = start, 
		end = end, thin = thin, verbose = verbose, n.report = n.report)
	theta <- rec$p.theta.recover.samples
	B <- rec$p.beta.recover.samples

	# extract other needed parameters
	sigmasq <- theta[,"sigma.sq"]
	phi <- theta[, "phi"]
	ev <- rep(0, nsim)
	fv <- rep(0, nsim)
	nu <- rep(0.5, nsim)

	# smoothness parameter for matern
	if(sp.obj$cov.model == "matern")
	{
		nu <- theta[, "nu"]
	}
	
	if(sp.obj$cov.model == "exponential")
	{
		cov.model <- 1
	}
	else if(sp.obj$cov.model == "gaussian")
	{
		cov.model <- 2
	}
	else if(sp.obj$cov.model == "matern")
	{
		cov.model <- 3
	}
	else
	{
		cov.model <- 4
	}
	
	# error.var and finescale.var parameters when nugget present
	if(colnames(theta)[2] == "tau.sq")
	{
		if(noisy)
		{
			fv <- theta[, "tau.sq"]		
		}else
		{
			ev <- theta[, "tau.sq"]		
		}
	}

	method.int <- 1
	if(method == "chol") method.int <- 2
	if(method == "svd") method.int <- 3
	
	out = .Call( "spLMPredict", ys = sp.obj$Y, 
		coordss = sp.obj$coords, pcoordss = pred.coords, 
		Xs = sp.obj$X, 
		Xps = pred.covars, 
		Bs = B, 
		sigmasqs = sigmasq, phis = phi, nus = nu, 
		evs = ev, fvs = fv, 
		cov_models = cov.model,
		methods = method.int, nreports = n.report, verboses = as.numeric(verbose), 
		PACKAGE = "SpatialTools")
  class(out) = "jointPredictiveSample"
  return(out)
}


# spLMPredictJoint <- function(sp.obj, pred.coords, pred.covars, start = 1, 
# 	end = nrow(sp.obj$p.theta.samples), thin = 1, verbose = TRUE, n.report = 100, 
# 	noisy = FALSE, method = "eigen")
# {
# 	nsim <- end - start + 1
# 	n <- nrow(sp.obj$coords)
# 	np <- nrow(pred.coords)
# 	opos <- 1:n
# 	ppos <- (n+1):(n+np)
# 
# 	rec <- spBayes::spRecover(sp.obj, get.beta = TRUE, get.w = FALSE, start = start, end = end, 
# 		thin = thin, verbose = verbose, n.report = n.report)
# 	theta <- rec$p.theta.recover.samples
# 	B <- rec$p.beta.recover.samples
# 
# 	# create other needed parameters
# 	error <- rep(0, nsim)
# 	finescale <- rep(0, nsim)
# 	nu <- rep(0.5, nsim)
# 
# 
# 	# smoothness parameter for matern
# 	if(sp.obj$cov.model == "matern")
# 	{
# 		nu <- theta[, "nu"]
# 	}
# 	
# 	# error.var and finescale.var parameters when nugget present
# 	if(colnames(theta)[2] == "tau.sq")
# 	{
# 		if(noisy)
# 		{
# 			finescale <- theta[, "tau.sq"]		
# 		}else
# 		{
# 			error <- theta[, "tau.sq"]		
# 		}
# 	}
# 	
# 	yp.sim <- matrix(0, nrow = np, ncol = nsim)
# 
# 	D <- dist1(sp.obj$coords)
# 	Dp <- dist1(pred.coords)
# 	Dop <- dist2(sp.obj$coords, pred.coords)
# 
# 	if(verbose)
# 	{
# 		cat("Samples from joint posterior: ")
# 	}
# 	for(i in 1:nsim)
# 	{
# 		mu <- sp.obj$X %*% B[i,]
# 		mup <- pred.covars %*% B[i,]
# 
# 		Va <- cov.sp(sp.obj$coords, sp.type = sp.obj$cov.model, 
# 			sp.par = c(theta[i,"sigma.sq"], 1/theta[i, "phi"]), 
# 			error.var = error[i], smoothness = nu[i], 
# 			finescale.var = finescale[i], pcoords = pred.coords, D = D, Dp = Dp, Dop = Dop)
# 
# 		yp.sim[, i] <- rcondnorm(1, y = sp.obj$Y, mu = mu, mup = mup, 
# 			V = Va$V, Vp = Va$Vp, Vop = Va$Vop, method = method)
# 		
# 		if(verbose)
# 		{	
# 			if(i %% n.report == 0){ cat(paste(i,""))}
# 		}
# 	}
# 	if(verbose)
# 	{
# 		cat(paste("\n"))
# 	}
# 	return(yp.sim)
# }
