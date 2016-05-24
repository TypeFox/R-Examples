msBP.Gibbs <-
function(x, a, b, g0 = "normal", g0par=c(0,1), mcmc, grid = list(n.points=40, low=0.001, upp=0.999), state=NULL, hyper=NULL, printing=0, maxScale=5, ...){
	n = length(x)
	x = sort(x)
	x.grid = seq(grid$low,grid$upp, length=grid$n.points)
	if(is.null(hyper$hyperprior$a)) 
	{
		hyper$hyperprior$a <- FALSE
		hyper$hyperpar$beta <- FALSE
		hyper$hyperpar$gamma <- FALSE
	}
	else
	{
		if((hyper$hyperprior$a==TRUE) & (is.null(hyper$hyperpar$beta)  | is.null(hyper$hyperpar$beta) ) ) 
			stop("If a is random, insert values for beta and gamma hyperparamters\n")
		if((hyper$hyperprior$a==FALSE)) 
		{
			hyper$hyperpar$beta <- FALSE
			hyper$hyperpar$gamma <- FALSE
		}
	}	
	if(is.null(hyper$hyperprior$b))
	{
		hyper$hyperprior$b <- FALSE
		hyper$hyperpar$delta <- FALSE
		hyper$hyperpar$lambda <- FALSE
		hyper$hyperpar$gridB <- FALSE
	}
	else
	{
		if((hyper$hyperprior$b==TRUE) & (is.null(hyper$hyperpar$delta)  | is.null(hyper$hyperpar$lambda) ) ) 
			stop("If b is random, insert values for delta and lambda hyperparamters\n")
		if((hyper$hyperprior$b==FALSE)) 
		{
			hyper$hyperpar$delta <- FALSE
			hyper$hyperpar$lambda <- FALSE
			hyper$hyperpar$gridB <- FALSE
		}
	}
	if(is.null(hyper$hyperprior$g0))
	{
		hyper$hyperprior$g0 <- FALSE
		hyper$hyperpar$mu0 <- FALSE
		hyper$hyperpar$kappa0 <- FALSE
		hyper$hyperpar$alpha0 <- FALSE
		hyper$hyperpar$beta0 <- FALSE		
	}
	else
	{
		if((hyper$hyperprior$g0==TRUE) & (is.null(hyper$hyperpar$mu0)  | is.null(hyper$hyperpar$kappa0) | is.null(hyper$hyperpar$alpha0) | is.null(hyper$hyperpar$beta0) ) ) 
			stop("If hyperprior for G0 is assumed, please insert values for mu0, kappa0, alpha0, and beta0\n")
			if((hyper$hyperprior$g0==FALSE)) 
			{
				hyper$hyperpar$mu0 <- FALSE
				hyper$hyperpar$kappa0 <- FALSE
				hyper$hyperpar$alpha0 <- FALSE
				hyper$hyperpar$beta0 <- FALSE	
			}
	}
	if((hyper$hyperprior$g0==TRUE) & (g0 != "normal") )  
		stop("Hyperprior for the base G0 can be espressed only for g0='normal'\n")
	
	if((hyper$hyperprior$g0==TRUE) & (g0 == "normal") )
		{
			y = pnorm(x, g0par[1], g0par[2])
			y.grid = x.grid
			g0_x = rep(1,grid$n.points)
			g0_xi = rep(1,n)
		}
	else
	{		  
		if(g0 == "normal")
		{
			y = pnorm(x, g0par[1], g0par[2])
			y.grid = pnorm(x.grid, g0par[1], g0par[2])
			g0_x = dnorm(x.grid, g0par[1], g0par[2])
			g0_xi = dnorm(x, g0par[1], g0par[2])
		}
		if(g0 == "unif")
		{
			y = x
			y.grid = x.grid
			g0_x = rep(1,grid$n.points)
			g0_xi = rep(1,n)
		}
		if(g0 == "gamma")
		{
			y = pgamma(x, g0par[1], g0par[2])
			y.grid = pgamma(x.grid, g0par[1], g0par[2])
			g0_x = dgamma(x.grid, g0par[1], g0par[2])
			g0_xi = dgamma(x, g0par[1], g0par[2])
		}
		if(g0 == "empirical")
		{
			grid=list(n.points=grid$n.points, low=grid$low, upp=grid$upp)
			x.grid = seq(grid$low,grid$upp, length=grid$n.points)
			kern.smooth <- density(x, from=grid$low, to=grid$upp, n=grid$n.points)
			mass <- sum(kern.smooth$y)*mean(diff(kern.smooth$x))
			y.grid = cumsum(kern.smooth$y)*mean(diff(kern.smooth$x))
			g0_x = kern.smooth$y
			y = approxfun(x.grid, y.grid)(x)	
			g0_xi = approxfun(x.grid, g0_x)(x)	
		}
		if((g0 != "unif") & (g0 != "normal") & (g0 != "gamma") & (g0 != "empirical"))
			stop("Only normal, uniform, gamma and empirical Bayes allowed for version 1.1\n")
	}
	if(is.null(state))
	{
		state = list()
		clusters <- msBP.leaf.allocation(y, maxScale)
		state$sclus = clusters$s
		state$hclus = clusters$h
		startTrees = msBP.rtree(a, b, maxScale)
		state$Rstart = startTrees$R
		state$Sstart = startTrees$S
		state$wstart = msBP.compute.prob(startTrees)
	}
	res <- .C("msBPgibbs", 
		x=as.double(x),
		y=as.double(y), 
		par=as.double(c(a,b,hyper$hyperpar$beta, hyper$hyperpar$gamma, 
			hyper$hyperpar$delta, hyper$hyperpar$lambda, 
			hyper$hyperpar$mu0, hyper$hyperpar$kappa0, hyper$hyperpar$alpha0, hyper$hyperpar$beta0)), 
		sclus = as.integer(state$sclus),
		hclus = as.integer(state$hclus),		
		Sstart = as.double(tree2vec(state$Sstart)),
		Rstart = as.double(tree2vec(state$Rstart)),
		wstart = as.double(tree2vec(state$wstart)),
		hyperpar=as.integer(c(hyper$hyperprior$a,hyper$hyperprior$b,hyper$hyperprior$g0)), 
		nrep=as.integer(mcmc$nrep), 
		nb=as.integer(mcmc$nb), 
		aux=as.integer(c(n, maxScale, (2^(maxScale+1)-1), state$wstart$max.s)),
		printing = as.integer(c(mcmc$ndisplay,printing)), 
		grid=as.double(y.grid), 
		ngrid=as.integer(grid$n.points),
		griddyB = as.double(hyper$hyperpar$gridB),
		griddy_length = as.integer(length(hyper$hyperpar$gridB)),
		postDens=as.double(rep(0,mcmc$nrep*grid$n.points)), 
		postScale=as.double(rep(0, (maxScale+1)*mcmc$nrep)), 
		postS=as.double(rep(0, mcmc$nrep*(2^(maxScale+1)-1))), 
		postR=as.double(rep(0, mcmc$nrep*(2^(maxScale+1)-1))), 
		postpi=as.double(rep(0, mcmc$nrep*(2^(maxScale+1)-1))), 
		postA=as.double(rep(a, mcmc$nrep)), 
		postB=as.double(rep(b, mcmc$nrep)), 
		posts=as.integer(rep(0, n*mcmc$nrep)),
		posth=as.integer(rep(0, n*mcmc$nrep)), 
		postmu=as.double(rep(g0par[1], mcmc$nrep)),
		postsigma2=as.double(rep((g0par[2]^2), mcmc$nrep)),
		PACKAGE = "msBP"
	)
	cat("Iteration", mcmc$nrep, "over", mcmc$nrep, "\n")
	postDens <- matrix(res$postDens, nrow=res$nrep, ncol=res$ngrid, byrow=TRUE)[(mcmc$nb+1):mcmc$nrep,]
	postDens <- t(t(postDens) * g0_x)	
	postMeanDens <- apply(postDens, 2, mean)
	postDensSort <- apply(postDens, 2, sort)
	postLowDens <- postDensSort[(mcmc$nrep-mcmc$nb)*0.025,]
	postUppDens <- postDensSort[(mcmc$nrep-mcmc$nb)*0.975,]
	scale <- matrix(res$postScale, nrow=res$nrep, ncol=maxScale+1, byrow=TRUE)[(mcmc$nb+1):mcmc$nrep,]
	postMeanScale <- apply(scale, 2, mean)
	postS <- matrix(res$postS, nrow=res$nrep, ncol=(2^(maxScale+1)-1), byrow=TRUE)[(mcmc$nb+1):mcmc$nrep,]
	postMeanS <- apply(postS, 2, mean)
	postMeanS <- vec2tree(postMeanS)
	postR <- matrix(res$postR, nrow=res$nrep, ncol=(2^(maxScale+1)-1), byrow=TRUE)[(mcmc$nb+1):mcmc$nrep,]
	postMeanR <- apply(postR, 2, mean)
	postMeanR <- vec2tree(postMeanR)
	postW <- matrix(res$postpi, nrow=res$nrep, ncol=(2^(maxScale+1)-1), byrow=TRUE)[(mcmc$nb+1):mcmc$nrep,]
	postMeanW <- apply(postW, 2, mean)
	postMeanW <- vec2tree(postMeanW)
	posts <- matrix(res$posts, nrow=res$nrep, ncol=n, byrow=TRUE)[(mcmc$nb+1):mcmc$nrep,]
	postMeans <- apply(posts, 2, mean)
	posth <- matrix(res$posth, nrow=res$nrep, ncol=n, byrow=TRUE)[(mcmc$nb+1):mcmc$nrep,]
	postMeanh <- apply(posth, 2, mean)
 
 	density=list(postMeanDens=postMeanDens, postLowDens=postLowDens, postUppDens=postUppDens, xDens = x.grid)
	
 	mcmc.out=list(dens=postDens, a=res$postA[(mcmc$nb+1):mcmc$nrep], b=res$postB[(mcmc$nb+1):mcmc$nrep], 
		scale=scale, S=postS, R=postR, weights=postW, s=posts, h = posth, 
		mu=res$postmu[(mcmc$nb+1):mcmc$nrep], sigma2=res$postsigma2[(mcmc$nb+1):mcmc$nrep])
	if(hyper$hyperprior$g0!=TRUE) mcmc.out = mcmc.out[-c(10,11)]
	
	postmean = list(a=mean(res$postA[(mcmc$nb+1):mcmc$nrep]), b=mean(res$postB[(mcmc$nb+1):mcmc$nrep]),
		S=postMeanS, R=postMeanR, weights=postMeanW, scales=postMeanScale, 
	 	mu=mean(res$postmu[(mcmc$nb+1):mcmc$nrep]), sigma2=mean(res$postsigma2[(mcmc$nb+1):mcmc$nrep]))
	if(hyper$hyperprior$g0!=TRUE) postmean = postmean[-c(7,8)]

 	lpml <- msBP.LPML(y, g0_xi, postW, hyper$hyperprior$g0, res$postmu, res$postsigma2, x)
	
	list(density=density, mcmc=mcmc.out, postmean=postmean, fit=lpml)
}
