# Requires rjags:
if(require('rjags')){


	library('runjags')
	runjags.options(inits.warning=FALSE, rng.warning=FALSE, blockignore.warning=FALSE)

	library('parallel')

	testnum <- 1

	themodel <- "
	model{

		for(i in 1:N){
			Y[i] ~ dnorm(true.y[i], precision)
			true.y[i] <- (m * X[i]) + c
		}
		m ~ dunif(-1000,1000)
		c ~ dunif(-1000,1000)
		precision ~ dexp(1)

		#data# N, X
	}"

	# Simulate the data
	set.seed(1)
	N <- 20
	X <- 1:N
	Y <- rnorm(length(X), 2*X + 1, 1)

	# Some initial values to use for 2 chains:

	initfun <- function(chain){

		# data is made available within this function when it
		# is evaluated for each simulation:
		stopifnot(length(data$X) == data$N)

		m <- c(-10,10)[chain]
		c <- c(10,-10)[chain]
		precision <- c(0.01,100)[chain]

		.RNG.seed <- chain
		.RNG.name <- c("base::Super-Duper",
		"base::Wichmann-Hill")[chain]

		return(list(m=m, c=c, precision=precision,
		.RNG.seed=.RNG.seed, .RNG.name=.RNG.name))
	}

	# A simple function that removes (over-writes with NA) one datapoint at a time:
	datafun <- function(s){
		simdata <- Y
		simdata[s] <- NA
		return(list(Y=simdata))
	}


	# Set up a cluster to use with the parLapply method:
	cat('Running study test number', testnum, '\n'); testnum <- testnum+1
	cl <- makeCluster(2)
	# Call the simulations over the snow cluster:
	results <- run.jags.study(simulations=4, model=themodel, datafunction=datafun,
	targets=list(Y=Y, m=2, c=1), n.chains=2, inits=initfun, cl=cl)


	m <- 'model{
	d[1] ~ dpois(mu)
	d[2] ~ dpois(mu)
	d[3] ~ dpois(mu)
	mu ~ dgamma(1,1)
	#monitor# mu
	#data# d
	}'

	##### Can't test any more than 2 spawned processes on winbuilder #####

	cat('Running study test number', testnum, '\n'); testnum <- testnum+1
	mu <- list(1,1)
	d <- c(5, 4, 7)
	jr <- run.jags(m, method='rjags', n.chains=2, inits=list(list(mu=1), list(mu=1)), silent.jags=TRUE)
	# Drop 1 (would create 3 clusters except we pass it cl):
	r <- drop.k(jr, dropvars='d', cl=cl)
	stopCluster(cl)

	# Drop k (use lapply so we don't create a cluster with 4 nodes):
	cat('Running study test number', testnum, '\n'); testnum <- testnum+1
	r <- drop.k(jr, dropvars='d', simulations=4, k=2, silent.jags=TRUE, parallel.method=lapply)
	
	cat("All study/drop-k checks passed\n")
	
}else{
	cat("Note: the rjags package is not installed, so the study/drop-k tests were skipped\n")
}



