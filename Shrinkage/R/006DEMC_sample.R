
# Corey M. Yanofsky, 17 February 2009

DEMC_sample <- function(logtarget, n, states, blockindex, current_logprob, temperature_schedule, gamma_schedule,
			jitter_magnitude){
# DEMC_SAMPLE Perform Differential Evolution-Markov chain Monte Carlo sampling
# (see References)
#
# Usage
#
# states = DEMC_sample(logtarget, n,  initial_states, blockindex,...
#       initial_logprob, temperature_schedule, gamma_schedule, jitter_magnitude)
#
# INPUTS
#
# logtarget - a function accepting a vector of length d and returning the log
# probability density for that vector. The density may be unnormalized.
# (DEMC may fail if the target density does not have support over R^d.)
#
# n - number of iterates to sample
#
# initial_states - a list of length m containing states of the Markov chains in
# the population.
#
# OPTIONAL INPUTS 
#
# blockindex - a list of length b, each of the b entries containing the indexes
# of state vector elements belonging to the b'th block. For DE-MC without blocks,
# set blockindex = list(1:d); this is also the default.
#
# initial_logprob - a vector of length m containing the log probability of the
# initial states; if this has been pre-calculated then this can save time,
# especially if this routine is being called in a loop as part of a hybrid
# sampling scheme.
#
# temperature_schedule - schedule for temperatures used in simulated annealing.
# After the schedule is complete, the temperature is set to 1. An empty schedule
# defaults to a constant temperature of 1.
#
# gamma_schedule - schedule for scale factor for vector differences; the default value
# an unchanging schedule of 2.38/sqrt(2d), which is optimal for d-dimensional
# Gaussian targets. For example, a schedule of 9 iterations at gamma = 2.38/sqrt(2.*d)
# followed by one iteration at 0.99 is specified by
# gamma = c(rep(2.38/sqrt(2.*d),9),0.99)
#
# jitter_magnitude - a vector of length d containing scale factors for the jitter; the default
# value is rep(1e-5,d)
#
# OUTPUT - a list of length m containing the sampled states.
#
# References
#
# Cajo J.F. Ter Braak, "A Markov Chain Monte Carlo version of the genetic algorithm
# Differential Evolution: easy Bayesian computing for real parameter spaces"
# Stat Comput (2006) 16:239-249

# initialize and set defaults

	m <- length(states)
	d <- length(states[[1]])
	for(i in 1:m) stopifnot(d == length(states[[i]]))

	if(missing(jitter_magnitude)) jitter_magnitude <- rep(1e-5,d)
	if(missing(gamma_schedule)) gamma_schedule <- 2.38/sqrt(2*d)
	if(missing(temperature_schedule)) temperature_schedule <- rep(1,n)
	if(missing(current_logprob)) current_logprob <- sapply(states, logtarget)
	if(missing(blockindex)) blockindex <- list(1:d)
	nblocks <- length(blockindex);

	# fill out temperature schedule
	if(length(temperature_schedule) < n) temperature_schedule[(length(temperature_schedule)+1):n] <- 1

	# create gamma index
	gamma_schedule = rep(gamma_schedule, length.out = n)

	saved_states = vector("list",n)

	DEMC_chain_sample <- function() 
	{
		# sample two states other than i uniformly without replacement
		n1 <- ceiling(runif(1)*(m-1))
		n2 <- ceiling(runif(1)*(m-2))
		if(n1 >= i)  
		{
		    n1 <- n1 + 1
		    if(n2 >= i) n2 <- n2 + 1
		    if(n2 >= n1) n2 <- n2 + 1
		}
		else 
		{
		    if(n2 >= n1) n2 <- n2 + 1
		    if(n2 >= i) n2 <- n2 + 1
		}

		# generate proposal via differential evolution rule
		difference_vector <- rep(0,d)
		difference_vector[blockindex[[b]]] <- gamma_schedule[iter]*(states[[n1]][blockindex[[b]]] - states[[n2]][blockindex[[b]]]) + jitter_magnitude[blockindex[[b]]]*rnorm(length(blockindex[[b]]))
		proposal <- states[[i]] + difference_vector

		proposal_logprob <- logtarget(proposal)

		# Metropolis accept/reject step
		if(temperature_schedule[iter]*log(runif(1)) < (proposal_logprob - current_logprob[i]))
		{
	   	 states[[i]] <<- proposal
	   	 current_logprob[i] <<- proposal_logprob
		}
	} # end function

	# sample new states n times...
	for(iter in 1:n){
	    # for each of the m chains...
	    for(i in 1:m){
		# for each block
		for(b in 1:nblocks) {
			DEMC_chain_sample()
		} # end blocks
	    } # end chains
	    saved_states[[iter]] <- states;
	}# end iterates
	saved_states
}# end function
