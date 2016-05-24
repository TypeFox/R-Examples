## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(cache=TRUE, autodep=TRUE, fig.width=5, fig.height=5)

## ----package_loading, hide=TRUE, error=FALSE, warning=FALSE, message=FALSE----
library(epiR)  # For the BetaBuster function
library(compiler)  # To compile the larger functions for computational speed
library(coda)  # For processing Bayesian model output
library(shape)  # For nice colorbar legends
library(scales)  # For transparent colors
library(EpiBayes)  # Load our package

## ----betabuster_test-----------------------------------------------------
epi.betabuster(0.10, 0.95, FALSE, 0.15)

## ----example1_run, eval=TRUE---------------------------------------------
set.seed(2015)  # To ensure reproducible results
example1.run = EpiBayes_s(
    H = 1,  # 1 subzone
    k = rep(40, 1),  # 40 farms total
    n = c(rep(100, 35), rep(500, 5)),  #100 
    	# mollusks sampled in 35 farms and 500 sampled in the remaining 5 farms
    seasons = rep(c(1, 2, 3, 4), each = 10),  
    	# Seasons corresponding to each cluster 
	# (1 for  summer,  2 for  fall,  3 for  winter,  4 for  spring)
    	mumodes = matrix(c(
		0.50, 0.70, 
		0.50, 0.70, 
		0.02, 0.50, 
		0.02, 0.50
		), 4, 2, byrow = TRUE
	),  # Modes and 95th percentiles of 
	#subject - level prevalences for  each season in order
    reps = 10,  # 10 replicated data sets in this simulation
    MCMCreps = 100,  # 100 MCMC iterations per replicated data set 
    	# (increasing this would be a good idea for  real data but 
	# slows things down a lot) 
    poi = "tau",  # Want inference on cluster-level prevalence
    y = NULL,  # Leave this as NULL if we are doing simulation and not posterior 
    	# inference with a particular data set
    pi.thresh = 0.10,  # Have a 10% within-cluster design prevalence
    tau.thresh = 0.02,  # Have a 2% between-cluster design prevalence
    gam.thresh = 0.10,  # Doesn't matter since we have a 2-level model
    tau.T = 0.20,  # The "true cluster - level prevalence" that we simulate our data 
    	# with (this means about 20% of our clusters in each replicated data set 
		# will be diseased and will have a truly positive subject - 
		# level prevalence)
    poi.lb = 0,  # The lower bound for  estimating the cluster - level 
    	# prevalence (not of interest here)
    poi.ub = 1,  # The upper bound for  estimating the cluster - level 
    	# prevalence (not of interest here)
    p1 = 0.95,  # The probability (used like a confidence) that we must show our 
    	# cluster - level prevalence is above 2% in order to count that replicated 
	# data set as one in which we detected the disease
    psi = 4,  # The variability of the prevalences among infected clusters 
    	# within the subzone
    omegaparm = c(1000, 1),  # Prior parameters for omegamat (the probability 
    	# of the disease being in the region)
    gamparm = c(1000, 1),  # Prior parameters for  gammat (the subzone-level 
    	# prevalence)
    tauparm = c(1, 1),  # Prior parameters for  taumat (the cluster - level 
    	# prevalence)
    etaparm = c(15, 130),  # Prior parameters for  etamat (the diagnostic test 
    	# sensitivity)
    thetaparm = c(1000, 1),  # Prior parameters for  thetamat (the diagnostic 
    	# test specificity)
    burnin = 10  # The amount of MCMC iterations to "burn"
    )


## ----example1_p4, eval = TRUE--------------------------------------------
## Summary
example1.sum = summary(example1.run, n.output = 5)
example1.sum

## Plot the posterior distributions of cluster-level prevalence 
plot(example1.run)


## ----example1_traceplots, eval = TRUE------------------------------------
## Tau
	# Tau for  first replicated data set
	plot(example1.run$taumat[1, 1, ], type = "l") 
	# Tau for  second replicated data set
	plot(example1.run$taumat[2, 1, ], type = "l") 
	# Tau for  third replicated data set
	plot(example1.run$taumat[3, 1, ], type = "l") 
	# Tau for  tenth replicated data set
	plot(example1.run$taumat[10, 1, ], type = "l") 

## Pi
	# Pi for  the first farm (100 mollusks) in the first replicated data set
	plot(example1.run$pimat[1, 1, ], type = "l") 
	# Pi for  the tenth farm (100 mollusks) in the first replicated data set
	plot(example1.run$pimat[1, 10, ], type = "l") 
	# Pi for  the thirty-eighth farm (500 mollusks) in the first replicated data set
	plot(example1.run$pimat[1, 38, ], type = "l") 
	# Pi for  the first farm (100 mollusks) in the tenth replicated data set
	plot(example1.run$pimat[10, 1, ], type = "l") 
	# Pi for  the tenth farm (100 mollusks) in the tenth replicated data set
	plot(example1.run$pimat[10, 10, ], type = "l")


## ----example1_histograms, eval=TRUE--------------------------------------
par(mfrow=c(2,2))
## Tau
	# Tau for  first replicated data set
	hist(example1.run$taumat[1, 1, ], col = "cyan");box("plot") 
	# Tau for  second replicated data set
	hist(example1.run$taumat[2, 1, ], col = "cyan");box("plot") 
	# Tau for  third replicated data set
	hist(example1.run$taumat[3, 1, ], col = "cyan");box("plot") 
	# Tau for  tenth replicated data set
	hist(example1.run$taumat[10, 1, ], col = "cyan");box("plot") 

## Pi
	# Pi for  the first farm (100 mollusks) in the first replicated data set
	hist(example1.run$pimat[1, 1, ], col = "cyan");box("plot") 
	# Pi for  the tenth farm (100 mollusks) in the first replicated data set
	hist(example1.run$pimat[1, 10, ], col = "cyan");box("plot") 
	# Pi for  the thirty-eighth farm (500 mollusks) in the first replicated data set
	hist(example1.run$pimat[1, 38, ], col = "cyan");box("plot") 
	# Pi for  the first farm (100 mollusks) in the tenth replicated data set
	hist(example1.run$pimat[10, 1, ], col = "cyan");box("plot") 
	# Pi for  the tenth farm (100 mollusks) in the tenth replicated data set
	hist(example1.run$pimat[10, 10, ], col = "cyan");box("plot") 


## ----example2_run, eval=TRUE---------------------------------------------
set.seed(2015) # To ensure reproducible results
example2.run = EpiBayes_ns(
   H = 1,  
    k = rep(40, 1),  
    n = c(rep(100, 35), rep(500, 5)),  
    seasons = rep(c(1, 2, 3, 4), each = 10),  
    mumodes = matrix(c(
	    0.50, 0.70, 
	    0.50, 0.70, 
	    0.02, 0.50, 
	    0.02, 0.50
	    ), 4, 2, byrow = TRUE
    ),
    reps = 1,  # Only 1 replication this time
    MCMCreps = 1000,  # 1000 MCMC iterations now that we have only 1 replication  
    poi = "tau",  # Want inference on cluster-level prevalence
    y = matrix(0, nrow = 1, ncol = 40),  # Set so we have all negative diagnostic 
    	# test results
    pi.thresh = 0.10,  # Have a 10% within-cluster design prevalence
    tau.thresh = 0.02,  # Doesn't matter since we aren't simulating
    gam.thresh = 0.10,  # Doesn't matter since we aren't simulating
    tau.T = 0.20,  # Doesn't matter since we aren't simulating any data in this 
    	# case (we'll leave it the same here)
    poi.lb = 0,  # Doesn't matter since we aren't simulating
    poi.ub = 1,  # Doesn't matter since we aren't simulating
    p1 = 0.95,  # Doesn't matter since we aren't simulating
    psi = 4,  
    omegaparm = c(1000, 1),  
    gamparm = c(1000, 1),  
    tauparm = c(1, 1),  
    etaparm = c(15, 130),  
    thetaparm = c(1000, 1),  
    burnin = 10
    )	


## ----example2_summary, eval = TRUE---------------------------------------
## Summary
example2.sum = summary(example2.run)
example2.sum

## Plot the posterior distributions of cluster-level prevalence 
plot(example2.run)


## ----example2_coda, eval=TRUE--------------------------------------------
## Must transform the posterior distribution of tau into an mcmc 
# object for  CODA to recognize it
example2.tauposterior = as.mcmc(example2.run$taumat[1, 1, ])

## Now we can use any function in the CODA package to analyze the posterior
	# Generic summary statistics 
	summary(example2.tauposterior)
	
	# 95% HPD interval for  the cluster - level prevalence
	HPDinterval(example2.tauposterior, prob = 0.95)
	
	# Geweke convergence statistic
	# If  it looks like the realization of a standard normal random variable
	coda::geweke.diag(example2.tauposterior)
	
	# Trace plot and density plot for  tau minus burnin
	plot(example2.tauposterior)

	# Autocorrelation plot for  the chain minus burnin
	autocorr.plot(example2.tauposterior)


## ----example3_run, eval = TRUE-------------------------------------------
example3.run = EpiBayesSampleSize(
	H = c(1, 2),  # Allow 1 and 2 subzones
	k = c(10, 20, 30),  # Allow 10, 20, and 30 clusters per subzone
	n = c(100, 300, 500),  # Allow 100, 300, and 500 subjects per 
		# cluster per subzone
	season = 3,  # Force all sampling to be done in Winter
	reps = 2, 
	MCMCreps = 2,
	tau.T = 0,
	y = NULL,
	poi = "tau",
	mumodes = matrix(c(
		0.50, 0.70, 
		0.50, 0.70, 
		0.02, 0.50, 
		0.02, 0.50
		), 4, 2, byrow = TRUE
	),
	pi.thresh = 0.10,
	tau.thresh = 0.02,
	gam.thresh = 0.10,
	poi.lb = 0,
	poi.ub = 0.1,
	p1 = 0.95,
	psi = 4,
	tauparm = c(1, 1),
	omegaparm = c(1000, 1),
	gamparm = c(1000, 1),
	etaparm = c(15, 130),
	thetaparm = c(100, 6),
	burnin = 1
	)


## ----example3_print, eval = TRUE-----------------------------------------
example3.run  # Prints all three simulation statistics

# Prints only those two requested
print(example3.run, out.ptilde = c("p4.tilde", "p6.tilde"))  


