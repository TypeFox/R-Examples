# test all of the runjags dispatch methods with a toy example:

library(runjags)

runjags.options(inits.warning=FALSE, rng.warning=FALSE, blockignore.warning=FALSE, silent.jags=TRUE)

# Require for as.mcmc.list and niter:
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


# Get the JAGS path - rjags adds the JAGS path to the PATH in Windows...
try(library(rjags))
jagspath <- findjags()

# Only run the JAGS methods tests if we have found JAGS and have permission to run it:
if(jagspath!="JAGS not found" && testjags(jagspath)$JAGS.available){
	
	testnum <- 1
	cat('Testing the simple method\n')
	# Try the simple method and if it doesn't work give a warning but don't fail (likely to be permissions problems)
	success <- try({
		results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='simple',temp=FALSE)
	})
	if(inherits(success, 'try-error')){
		cat("JAGS was found but the simple method failed; it is possible that there were permissions issues or similar.  Details as follows:\n")
		t <- testjags(silent=FALSE)
		cat(failed.jags('output')[[1]])
		print(file.info(jagspath))
		print(file.info(getwd())[,1:3])
		cat("All test methods (except possibly rjags) were skipped\n")
	}else{
		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		stopifnot(niter(as.mcmc.list(results))==1000)

		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='parallel')
		stopifnot(niter(as.mcmc.list(results))==1000)

		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='interruptible')
		stopifnot(niter(as.mcmc.list(results))==1000)
  
		# Same as in checkinputs but it's for rjags there:
		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results2 <- extend.jags(results, sample=1000, add.monitor="true.y", summarise=FALSE)
		stopifnot(nvar(as.mcmc.list(results2))==(3+N))
		stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[103]=='true.y[100]')

		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results2 <- extend.jags(results, sample=1000, add.monitor=c("true.y", "dic"), summarise=FALSE)
		stopifnot(nvar(as.mcmc.list(results2))==(3+N))
		stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[103]=='true.y[100]')

		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results2 <- extend.jags(results, sample=1000, add.monitor=c("true.y", "deviance"), summarise=FALSE)
		stopifnot(nvar(as.mcmc.list(results2))==(4+N))
		stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[104]=='deviance')
		
		# Same as in checkinputs but it's for rjags there:
		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results2 <- extend.jags(results, sample=1000, add.monitor="true.y", summarise=FALSE, method='simple')
		stopifnot(nvar(as.mcmc.list(results2))==(3+N))
		stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[103]=='true.y[100]')

		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results2 <- extend.jags(results, sample=1000, add.monitor=c("true.y", "dic"), summarise=FALSE)
		stopifnot(nvar(as.mcmc.list(results2))==(3+N))
		stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[103]=='true.y[100]')

		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results2 <- extend.jags(results, sample=1000, add.monitor=c("true.y", "deviance"), summarise=FALSE)
		stopifnot(nvar(as.mcmc.list(results2))==(4+N))
		stopifnot(varnames(as.mcmc.list(results2))[1]=='m' && varnames(as.mcmc.list(results2))[104]=='deviance')
		
		# Snow gives problems here ... but it does work!
		#results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='snow')
		#stopifnot(niter(as.mcmc.list(results))==1000)

		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		info <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='background')
		t <- 0
		repeat{
			# Change thin, chain and variables:
			s <- try(results <- results.jags(info, read.monitor='m', return.samples=100, recover.chains=2, summarise=FALSE))
			if(class(s)!='try-error') break
			if(t==30) stop("Timed out waiting for the bgparallel method")
			Sys.sleep(1)
			t <- t+1
		}
		stopifnot(niter(as.mcmc.list(results))==100)
		stopifnot(nvar(as.mcmc.list(results))==1)
		stopifnot(nchain(as.mcmc.list(results))==1)
	
	
		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		info <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='bgparallel')
		t <- 0
		repeat{
			s <- try(results <- results.jags(info))
			if(class(s)!='try-error') break
			if(t==30) stop("Timed out waiting for the bgparallel method")
			Sys.sleep(1)
			t <- t+1
		}
		stopifnot(niter(as.mcmc.list(results))==1000)	
				
		# Check combine.mcmc does what it says on the tin:
		stopifnot(niter(combine.mcmc(results, return.samples=1000, collapse.chains=TRUE))==1000)
		stopifnot(niter(combine.mcmc(results, return.samples=11, collapse.chains=TRUE))==11)
		stopifnot(niter(combine.mcmc(results, return.samples=100, collapse.chains=FALSE))==100)
		stopifnot(niter(combine.mcmc(results, thin=1, collapse.chains=TRUE))==2000)
		stopifnot(niter(combine.mcmc(results, thin=10, collapse.chains=TRUE))==200)
		stopifnot((niter(combine.mcmc(results, thin=15, collapse.chains=TRUE))*15)>=2000)
		stopifnot(niter(combine.mcmc(results, thin=10, collapse.chains=FALSE))==100)
		
		# Check we can use the extend wrapper:
		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		newres <- extend.jags(results, sample=0)
		stopifnot(newres$burnin==results$burnin)
		stopifnot(niter(as.mcmc.list(newres))==niter(as.mcmc.list(results)))
		
		# Check a single iteration works:
		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results <- run.jags(model, n.chains=2, sample=1, burnin=1000, inits=initfunction, method='interruptible', summarise=FALSE)
		stopifnot(niter(as.mcmc.list(results))==1)
		
		# And that precision can be ignored:
		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		results <- run.jags(model, monitor=c('m','c'), n.chains=2, sample=10, burnin=100, inits=initfunction, method='interruptible', summarise=FALSE, noread.monitor='precision')
		stopifnot(nvar(as.mcmc.list(results))==2)
		
		# Check the version number is correct:
		cat('Running method test number', testnum, '\n'); testnum <- testnum+1
		stopifnot(newres$runjags.version[1]==runjags:::runjagsprivate$runjagsversion)
		
	}
}else{
	cat("JAGS could not be called externally at the path: ", jagspath, "\n")
	cat("All test methods except possibly rjags and rjagsparallel were skipped\n")	
}

testnum <- 1
if(require("rjags")){
	cat('Running rjags method test number', testnum, '\n'); testnum <- testnum+1	
	results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='rjags')
	stopifnot(niter(as.mcmc.list(results))==1000)
	
	cat('Running rjags method test number', testnum, '\n'); testnum <- testnum+1
	results <- run.jags(model, n.chains=2, sample=1, burnin=1000, inits=initfunction, method='rjags', summarise=FALSE)
	stopifnot(niter(as.mcmc.list(results))==1)
	
	cat('Running rjags method test number', testnum, '\n'); testnum <- testnum+1
	runjags.options(silent.jags=TRUE, silent.runjags=TRUE, debug=FALSE)
	output <- capture.output(results <- run.jags(model, n.chains=2, sample=1000, burnin=1000, inits=initfunction, method='rjparallel'))
	stopifnot(length(output)==0)
	stopifnot(niter(as.mcmc.list(results))==1000)
}else{
	cat("The rjags and rjagsparallel method checks were skipped as the rjags package is not installed\n")
}

cat("All methods checks passed\n")

