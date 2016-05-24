# First load runjags and verify that the options can be set:
library('runjags')
runjags.options(inits.warning=FALSE, rng.warning=FALSE, silent.jags=TRUE)
newopts <- runjags.options()
stopifnot(newopts$inits.warning==FALSE)
stopifnot(newopts$rng.warning==FALSE)
stopifnot(newopts$silent.jags==TRUE)


# Test some runjags inputs with the rjags method:

testnum <- 1
if(require("rjags")){

	# Required for nvar etc:
	library("coda")

	model <- "model {
	for(i in 1 : N){ #data# N
	Y[i] ~ dnorm(true.y[i], precision); #data# Y
	true.y[i] <- (m * X[i]) + c; #data# X
	}
	m ~ dunif(-1000,1000); #inits# m
	c ~ dunif(-1000,1000);
	precision ~ dexp(1);
	#monitor# m, c, precision
	}"

	# Simulate the data
	X <- 1:100
	Y <- rnorm(length(X), 2*X + 10, 1)
	N <- length(X)

	initfunction <- function(chain) return(switch(chain, "1"=list(m=-10), "2"=list(m=10)))
	initfunction2 <- function() return(switch(sample(c(1,2),1), "1"=list(m=-10), "2"=list(m=10)))

	datalist <- list(X=X, Y=Y, N=N)
	
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='rjags')
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction2, method='rjags')
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=lapply(1:2,initfunction), method='rjags')
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction(1), method='rjags')

	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, data=datalist, inits=initfunction, method='rjags')
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, data=datalist, inits=initfunction2, method='rjags')
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, data=datalist, inits=lapply(1:2,initfunction), method='rjags')
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, data=datalist, inits=initfunction(1), method='rjags')
	stopifnot(nvar(as.mcmc.list(results))==(3))
	stopifnot(all(varnames(as.mcmc.list(results))==c('m','c','precision')))

	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- extend.jags(results,sample=1000)
	stopifnot(niter(as.mcmc.list(results))==2000)
	stopifnot(nvar(as.mcmc.list(results))==3)

	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results2 <- extend.jags(results, sample=1000, drop.chain=1, summarise=FALSE)
	stopifnot(nchain(as.mcmc.list(results2))==1)
	stopifnot(identical(list.format(results$end.state[[2]])$.RNG.name, list.format(results2$end.state[[1]])$.RNG.name))
  
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results2 <- extend.jags(results, sample=1000, drop.monitor="precision", summarise=FALSE)
	stopifnot(nvar(as.mcmc.list(results2))==2)
  
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results2 <- extend.jags(results, sample=1000, add.monitor="true.y", summarise=FALSE)
	stopifnot(nvar(as.mcmc.list(results2))==(3+N))
	stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[103]=='true.y[100]')

	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results2 <- extend.jags(results, sample=1000, add.monitor=c("true.y", "dic"), summarise=FALSE)
	stopifnot(nvar(as.mcmc.list(results2))==(3+N))
	stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[103]=='true.y[100]')

	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results2 <- extend.jags(results, sample=1000, add.monitor=c("true.y", "deviance"), summarise=FALSE)
	stopifnot(nvar(as.mcmc.list(results2))==(4+N))
	stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[104]=='deviance')
	
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, data=datalist, inits=initfunction, method='rjparallel')
	stopifnot(nvar(as.mcmc.list(results))==(3))
	stopifnot(all(varnames(as.mcmc.list(results))==c('m','c','precision')))	

	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results2 <- extend.jags(results, sample=1000, add.monitor="true.y", summarise=FALSE, method='rjp')
	stopifnot(nvar(as.mcmc.list(results2))==(3+N))
	stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[103]=='true.y[100]')
	
	cat('Running input test number', testnum, '\n'); testnum <- testnum+1
	results2 <- extend.jags(results, sample=1000, add.monitor=c("true.y", "deviance"), summarise=FALSE, method='rjp')
	stopifnot(nvar(as.mcmc.list(results2))==(4+N))
	stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[104]=='deviance')
	
	cat("All input checks passed\n")

}else{
	cat("The input checks were not performed as rjags is not installed\n")
}