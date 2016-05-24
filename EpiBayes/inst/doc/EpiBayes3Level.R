## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(cache=TRUE, autodep=TRUE, fig.width=5, fig.height=5)

## ----example1_obsmat, eval = TRUE----------------------------------------
obs.y = matrix(c(
		rep(0, 10),  # Subzone 1
		rep(0, 10),  # Subzone 2
		rep(10, 10), rep(15, 15), rep(0, 25)  # Subzone 3
		),
		nrow = 1
	)


## ----example1_run, eval=TRUE---------------------------------------------
set.seed(2015)  # To ensure reproducible results
example1.run = EpiBayes_s(
    H = 3,  # 3 subzones
    k = c(rep(10, 2), rep(50, 1)),  # 10 farms in two subzones, 50 in 
    	# the third subzone
    n = rep(100, 70),  #100 cows sampled in each of the 70 clusters
    seasons = rep(2, 70),  # Seasons corresponding to each cluster 
    	# (1 for  summer,  2 for  fall,  3 for  winter,  4 for  spring)
	mumodes = matrix(c(
		0.10, 0.50, 
		0.10, 0.50, 
		0.10, 0.50, 
		0.10, 0.50
		), 4, 2, byrow = TRUE
	), # Modes and 95th percentiles of 
    	# subject - level prevalences for  each season in order
    reps = 1,  # 1 replicated data set in this simulation
    MCMCreps = 100,  # 100 MCMC iterations per replicated data 
    	# set (increasing this would be a good idea for  real data but slows
	#  things down a lot) 
    poi = "tau",  # Want inference on cluster-level prevalence
    y = obs.y,  # Specify the number of positive test results we saw for each farm
    pi.thresh = 0.05,  # The 5% threshold (design prevalence) for  the 
    	# cluster - level prevalence 
    tau.thresh = 0.02,  # The 2% threshold (design prevalence) for  the 
    	# cluster - level prevalence 
    gam.thresh = 0.01,  # The 1% threshold (design prevalence) for  the 
    	# cluster - level prevalence 
    tau.T = 0.20,  # The "true cluster - level prevalence" that we simulate our 
    	# data with (this means about 20% of our clusters in each replicated 
	# data set will be diseased and will have a truly positive 
	# subject - level prevalence)
    poi.lb = 0,  # The lower bound for  estimating the cluster - level 
    	# prevalence (not of interest here)
    poi.ub = 1,  # The upper bound for  estimating the cluster - level 
    	# prevalence (not of interest here)
    p1 = 0.95,  # The probability (used like a confidence) that we must show 
    	# our cluster - level prevalence is above 2% in order to count that 
	# replicated data set as one in which we detected the disease
    psi = 4,  # The variability of the prevalences among infected clusters within 
    	# the subzone
    omegaparm = c(1000, 1),  # Prior parameters for omegamat (the probability 
    	# of the disease being in the region)
    gamparm = c(20, 30),  # Prior parameters for  gammat (the subzone-level 
    	# prevalence)
    tauparm = c(1, 1),  # Prior parameters for  taumat (the cluster - level 
    	# prevalence)
    etaparm = c(10, 1),  # Prior parameters for  etamat (the diagnostic 
    	# test sensitivity)
    thetaparm = c(10, 1),  # Prior parameters for  thetamat (the diagnostic 
    	# test specificity)
    burnin = 10  # The amount of MCMC iterations to "burn"
    )


## ----example1_summary, eval = TRUE---------------------------------------
## Summary
example1.sum = summary(example1.run)
example1.sum

## Plot the posterior distributions of cluster-level prevalence 
plot(example1.run)


## ----example1_traceplots, eval = TRUE------------------------------------
## Trace plots

## Tau
	# Tau for the first subzone
	plot(example1.run$taumat[1, 1, -c(1:10)], type = "l")
	# Tau for the second subzone
	plot(example1.run$taumat[1, 2, -c(1:10)], type = "l") 
	# Tau for the third subzone
	plot(example1.run$taumat[1, 3, -c(1:10)], type = "l") 

## Pi
	# Pi for the first farm in the first subzone
	plot(example1.run$pimat[1, 1, -c(1:10)], type = "l") 
	# Pi for the tenth farm in the first subzone
	plot(example1.run$pimat[1, 10, -c(1:10)], type = "l") 
	# Pi for the first farm in the second subzone
	plot(example1.run$pimat[1, 11, -c(1:10)], type = "l") 
	# Pi for the first farm in the third subzone
	plot(example1.run$pimat[1, 21, -c(1:10)], type = "l") 
	# Pi for the fiftieth farm in the third subzone
	plot(example1.run$pimat[1, 70, -c(1:10)], type = "l") 

## Histograms

## Tau
	# Tau for the first subzone
	plot(density(example1.run$taumat[1, 1, c(1:10)], from = 0, to = 1)) 

## Pi
	# Pi for the first farm in the first subzone
	plot(density(example1.run$pimat[1, 1, c(1:10)], from = 0, to = 1))


## ----example2_inputmat, eval = TRUE--------------------------------------
year = rep(c(2010:2014), each = 70)
subz = rep(rep(c("First", "Second", "Third"), c(10, 10, 50)), 5)
size = rep(100, 70*5)
season = rep(2, 70*5)
y = matrix(c(
		rep(0, 10), #Year 2010: Subzone 1
		rep(0, 10), #Year 2010: Subzone 2
		rep(10, 10), rep(15, 15), rep(0, 25), #Year 2010: Subzone 3
		rep(2, 10), #Year 2011: Subzone 1
		rep(0, 10), #Year 2011: Subzone 2
		rep(5, 10), rep(10, 15), rep(0, 25), #Year 2011: Subzone 3
		rep(0, 10), #Year 2012: Subzone 1
		rep(4, 10), #Year 2012: Subzone 2
		rep(0, 10), rep(5, 15), rep(0, 25), #Year 2012: Subzone 3
		rep(8, 10), #Year 2013: Subzone 1
		rep(0, 10), #Year 2013: Subzone 2
		rep(0, 10), rep(0, 15), rep(0, 25), #Year 2013: Subzone 3
		rep(4, 10), #Year 2014: Subzone 1
		rep(0, 10), #Year 2014: Subzone 2
		rep(0, 10), rep(0, 15), rep(0, 25) #Year 2014: Subzone 3
		),
		ncol = 1
	)
	
example2.inputdf = data.frame(year, subz, size, season, y)


## ----example2_run, eval=TRUE---------------------------------------------
set.seed(2015)
example2.run = EpiBayesHistorical(
	input.df = example2.inputdf,  # Our input matrix
	orig.tauparm = c(1, 1),  # tau prior parameters in the first year
	burnin = 1,  # Number of MCMC iterations to burn
	MCMCreps = 10,  # Number of MCMC iterations
	tau.T = 0.2,  # Doesn't matter since reps = 1
	poi = "tau",  # Leave parameter of interest as cluster-level prevalence
	mumodes = matrix(c(
		 0.10, 0.50, 
		 0.10, 0.50, 
		 0.10, 0.50, 
		 0.10, 0.50
		 ), 4, 2, byrow = TRUE
	 ),# Season-specific average subject-level 
		# prevalences in infected clusters 
         pi.thresh = 0.05,  # The 5% threshold (design prevalence) for  the 
    		# cluster - level prevalence 
         tau.thresh = 0.02,  # Doesn't matter since reps = 1
         gam.thresh = 0.01,  # Doesn't matter since reps = 1
	poi.lb = 0,  # Doesn't matter since reps = 1
	poi.ub = 1,  # Doesn't matter since reps = 1
	p1 = 0.95,  # Doesn't matter since reps = 1
	psi = 4,  # (related to) variability of subject-level prevalences in 
		# infected clusters
	omegaparm = c(1000, 1),  # Prior parameters for the probability of the
		# disease being in the region (almost always 1)
	gamparm = c(20, 30),  # Prior parameters for the subzone-level prevalence 
		# (mean of about 0.4)
	etaparm = c(10, 1),  # Prior parameters for diagnostic test sensitivity 
		# (mean of about 0.9)
	thetaparm = c(10, 1)  # Prior parameters for diagnostic test specificity 
		# (mean of about 0.9)
	)
	

## ----example2_plot, eval = TRUE------------------------------------------
## Plot the posterior distributions of cluster-level prevalence with 
	# one plotting window for each subzone and colors ranging 
	# across years
plot(example2.run)


## ----example2_summary, eval=TRUE-----------------------------------------
## Summaries
	# By mean
	example2.meansum = summary(example2.run, sumstat = "mean", 
		time.labels = 2010:2014)
	example2.meansum
	
	# By 95th percentiles
	## Summaries
	example2.95persum = summary(example2.run, sumstat = "quantile", 
		prob = 0.95, time.labels = 2010:2014)
	example2.95persum
	
## Plotting the summaries across time
	# Plot means
	plot(example2.meansum)
	# Can add a line to compare to a certain design prevalence
	abline(h = 0.05, lty = 2, col = "black", lwd = 2)  
	
	# Plot 95th percentiles
	plot(example2.95persum)
	# Can add a line to compare to a certain design prevalence
	abline(h = 0.05, lty = 2, col = "black", lwd = 2)  


