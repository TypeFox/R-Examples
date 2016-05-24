generateControlFile <- function(file = "controlfile.txt", type = "diversification", params = NULL) {
	
	templates <- list(diversification = 
	'# BAMM configuration file for speciation/extinction analysis 
	# ==========================================================
	#
	# Format
	# ------
	#
	#     - Each option is specified as: option_name = option_value
	#     - Comments start with # and go to the end of the line
	#     - True is specified with "1" and False with "0"
	
	
	################################################################################
	# GENERAL SETUP AND DATA INPUT
	################################################################################
	
	modeltype = speciationextinction        
	# Specify "speciationextinction" or "trait" analysis
	                                  
	treefile = %%%%
	# File name of the phylogenetic tree to be analyzed
	
	runInfoFilename = run_info.txt
	# File name to output general information about this run
	
	sampleFromPriorOnly = 0                 
	# Whether to perform analysis sampling from prior only (no likelihoods computed)
	
	runMCMC = 1                             
	# Whether to perform the MCMC simulation. If runMCMC = 0, the program will only
	# check whether the data file can be read and the initial likelihood computed
	
	loadEventData = 0                       
	# Whether to load a previous event data file
	
	eventDataInfile = event_data_in.txt
	# File name of the event data file to load, used only if loadEventData = 1
	
	initializeModel = 1                     
	# Whether to initialize (but not run) the MCMC. If initializeModel = 0, the
	# program will only ensure that the data files (e.g., treefile) can be read
	
	useGlobalSamplingProbability = 1        
	# Whether to use a "global" sampling probability. If False (0), expects a file
	# name for species-specific sampling probabilities (see sampleProbsFilename)
	                                        
	globalSamplingFraction = 1.0            
	# The sampling probability. If useGlobalSamplingFraction = 0, this is ignored
	# and BAMM looks for a file name with species-specific sampling fractions
	
	sampleProbsFilename = sample_probs.txt
	# File name containing species-specific sampling fractions
	
	# seed = 12345
	# Seed for the random number generator. 
	# If not specified (or is -1), a seed is obtained from the system clock
	
	overwrite = 0
	# If True (1), the program will overwrite any output files in the current
	# directory (if present)
	
	
	################################################################################
	# PRIORS
	################################################################################
	
	expectedNumberOfShifts = 1.0
	# prior on the number of shifts in diversification
	# Suggested values: 
	#     expectedNumberOfShifts = 1.0 for small trees (< 500 tips)
	#	  expectedNumberOfShifts = 10 or even 50 for large trees (> 5000 tips) 
	 
	lambdaInitPrior = 1.0
	# Prior (rate parameter of exponential) on the initial lambda value for rate
	# regimes
	
	lambdaShiftPrior = 0.05
	# Prior (std dev of normal) on lambda shift parameter for rate regimes
	# You cannot adjust the mean of this distribution (fixed at zero, which is
	# equal to a constant rate diversification process)
	
	muInitPrior = 1.0
	# Prior (rate parameter of exponential) on extinction rates  
	
	lambdaIsTimeVariablePrior = 1
	# Prior (probability) of the time mode being time-variable (vs. time-constant)
	
	
	################################################################################
	# MCMC SIMULATION SETTINGS & OUTPUT OPTIONS
	################################################################################
	
	numberOfGenerations = %%%%
	# Number of generations to perform MCMC simulation
	
	mcmcOutfile = mcmc_out.txt
	# File name for the MCMC output, which only includes summary information about
	# MCMC simulation (e.g., log-likelihoods, log-prior, number of processes)
	
	mcmcWriteFreq = 1000
	# Frequency in which to write the MCMC output to a file
	
	eventDataOutfile = event_data.txt
	# The raw event data (these are the main results). ALL of the results are
	# contained in this file, and all branch-specific speciation rates, shift
	# positions, marginal distributions etc can be reconstructed from this output.
	# See R package BAMMtools for working with this output
	
	eventDataWriteFreq = 1000
	# Frequency in which to write the event data to a file
	
	printFreq = 1000
	# Frequency in which to print MCMC status to the screen
	
	acceptanceResetFreq = 1000
	# Frequency in which to reset the acceptance rate calculation
	# The acceptance rate is output to both the MCMC data file and the screen
	
	# outName = BAMM
	# Optional name that will be prefixed on all output files (separated with "_")
	# If commented out, no prefix will be used
	
	
	################################################################################
	# OPERATORS: MCMC SCALING OPERATORS
	################################################################################
	
	updateLambdaInitScale = 2.0
	# Scale parameter for updating the initial speciation rate for each process
	
	updateLambdaShiftScale = 0.1
	# Scale parameter for the exponential change parameter for speciation
	
	updateMuInitScale = 2.0
	# Scale parameter for updating initial extinction rate for each process
	
	updateEventLocationScale = 0.05
	# Scale parameter for updating LOCAL moves of events on the tree
	# This defines the width of the sliding window proposal
	 
	updateEventRateScale = 4.0
	# Scale parameter (proportional shrinking/expanding) for updating
	# the rate parameter of the Poisson process 
	
	
	################################################################################
	# OPERATORS: MCMC MOVE FREQUENCIES
	################################################################################
	
	updateRateEventNumber = 0.1
	# Relative frequency of MCMC moves that change the number of events
	
	updateRateEventPosition = 1
	# Relative frequency of MCMC moves that change the location of an event on the
	# tree
	
	updateRateEventRate = 1
	# Relative frequency of MCMC moves that change the rate at which events occur 
	
	updateRateLambda0 = 1
	# Relative frequency of MCMC moves that change the initial speciation rate
	# associated with an event
	
	updateRateLambdaShift = 1
	# Relative frequency of MCMC moves that change the exponential shift parameter
	# of the speciation rate associated with an event
	
	updateRateMu0 = 1
	# Relative frequency of MCMC moves that change the extinction rate for a given
	# event
	
	updateRateLambdaTimeMode = 0
	# Relative frequency of MCMC moves that flip the time mode
	# (time-constant <=> time-variable)
	
	localGlobalMoveRatio = 10.0
	# Ratio of local to global moves of events 
	
	
	################################################################################
	# INITIAL PARAMETER VALUES
	################################################################################
	
	lambdaInit0 = 0.032
	# Initial speciation rate (at the root of the tree)
	
	lambdaShift0 = 0
	# Initial shift parameter for the root process
	
	muInit0 = 0.005
	# Initial value of extinction (at the root)
	
	initialNumberEvents = 0
	# Initial number of non-root processes
	 
	
	################################################################################
	# METROPOLIS COUPLED MCMC
	################################################################################
	
	numberOfChains = 4
	# Number of Markov chains to run
	
	deltaT = 0.1
	# Temperature increment parameter. This value should be > 0
	# The temperature for the i-th chain is computed as 1 / [1 + deltaT * (i - 1)]
	
	swapPeriod = 1000
	# Number of generations in which to propose a chain swap
	
	chainSwapFileName = chain_swap.txt
	# File name in which to output data about each chain swap proposal.
	# The format of each line is [generation],[rank_1],[rank_2],[swap_accepted]
	# where [generation] is the generation in which the swap proposal was made,
	# [rank_1] and [rank_2] are the chains that were chosen, and [swap_accepted] is
	# whether the swap was made. The cold chain has a rank of 1.
	
	
	################################################################################
	# NUMERICAL AND OTHER PARAMETERS
	################################################################################
	
	minCladeSizeForShift = 1
	# Allows you to constrain location of possible rate-change events to occur
	# only on branches with at least this many descendant tips. A value of 1
	# allows shifts to occur on all branches. 
	
	segLength = 0.02
	# Controls the "grain" of the likelihood calculations. Approximates the
	# continuous-time change in diversification rates by breaking each branch into
	# a constant-rate diversification segments, with each segment given a length
	# determined by segLength. segLength is in units of the root-to-tip distance of
	# the tree. So, if the segLength parameter is 0.01, and the crown age of your
	# tree is 50, the "step size" of the constant rate approximation will be 0.5.
	# If the value is greater than the branch length (e.g., you have a branch of
	# length < 0.5 in the preceding example) BAMM will not break the branch into
	# segments but use the mean rate across the entire branch.', 
	trait = '
	# BAMM configuration file for phenotypic analysis
	# ===============================================
	#
	# Format
	# ------
	#
	#     - Each option is specified as: option_name = option_value
	#     - Comments start with # and go to the end of the line
	#     - True is specified with "1" and False with "0"
	
	
	################################################################################
	# GENERAL SETUP AND DATA INPUT
	################################################################################
	
	modeltype = trait        
	# Specify "speciationextinction" or "trait" analysis
	                                  
	treefile = %%%%
	# File name of the phylogenetic tree to be analyzed
	
	traitfile = %%%%
	# File name of the phenotypic traits file
	
	runInfoFilename = run_info.txt
	# File name to output general information about this run
	
	sampleFromPriorOnly = 0                 
	# Whether to perform analysis sampling from prior only (no likelihoods computed)
	
	runMCMC = 1                             
	# Whether to perform the MCMC simulation. If runMCMC = 0, the program will only
	# check whether the data file can be read and the initial likelihood computed
	
	loadEventData = 0                       
	# Whether to load a previous event data file
	
	eventDataInfile = event_data_in.txt
	# File name of the event data file to load, used only if loadEventData = 1
	
	initializeModel = 1                     
	# Whether to initialize (but not run) the MCMC. If initializeModel = 0, the
	# program will only ensure that the data files (e.g., treefile) can be read
	
	# seed = 12345
	# Seed for the random number generator. 
	# If not specified (or is -1), a seed is obtained from the system clock
	
	overwrite = 0
	# If True (1), the program will overwrite any output files in the current
	# directory (if present)
	
	
	################################################################################
	# PRIORS
	################################################################################
	
	expectedNumberOfShifts = 1.0
	# prior on the number of shifts in diversification
	# Suggested values: 
	#     expectedNumberOfShifts = 1.0 for small trees (< 500 tips)
	#	  expectedNumberOfShifts = 10 or even 50 for large trees (> 5000 tips) 
	 
	betaInitPrior = 1.0
	# Prior (rate parameter of exponential) on the initial
	# phenotypic evolutionary rate associated with regimes
	
	betaShiftPrior = 0.05
	# Prior (std dev of normal) on the rate-change parameter
	# You cannot adjust the mean of this distribution (fixed at zero, which is
	# equal to a constant rate diversification process)
	
	useObservedMinMaxAsTraitPriors = 1
	# If True (1), will put a uniform prior density on the distribution
	# of ancestral character states, with upper and lower bounds determined
	# by the min and max of the observed data
	
	traitPriorMin = 0
	# User-defined minimum value for the uniform density on the distribution of
	# ancestral character states. Only used if useObservedMinMaxAsTraitPriors = 0.
	
	traitPriorMax = 0
	# User-defined maximum value for the uniform density on the distribution of
	# ancestral character states. Only used if useObservedMinMaxAsTraitPriors = 0.
	
	betaIsTimeVariablePrior = 1
	# Prior (probability) of the time mode being time-variable (vs. time-constant)
	
	
	################################################################################
	# MCMC SIMULATION SETTINGS & OUTPUT OPTIONS
	################################################################################
	
	numberOfGenerations = %%%%
	# Number of generations to perform MCMC simulation
	
	mcmcOutfile = mcmc_out.txt
	# File name for the MCMC output, which only includes summary information about
	# MCMC simulation (e.g., log-likelihoods, log-prior, number of processes)
	
	mcmcWriteFreq = 1000
	# Frequency in which to write the MCMC output to a file
	
	eventDataOutfile = event_data.txt
	# The raw event data (these are the main results). ALL of the results are
	# contained in this file, and all branch-specific speciation rates, shift
	# positions, marginal distributions etc can be reconstructed from this output.
	# See R package BAMMtools for working with this output
	
	eventDataWriteFreq = 1000
	# Frequency in which to write the event data to a file
	
	printFreq = 1000
	# Frequency in which to print MCMC status to the screen
	
	acceptanceResetFreq = 1000
	# Frequency in which to reset the acceptance rate calculation
	# The acceptance rate is output to both the MCMC data file and the screen
	
	# outName = BAMM
	# Optional name that will be prefixed on all output files (separated with "_")
	# If commented out, no prefix will be used
	
	
	################################################################################
	# OPERATORS: MCMC SCALING OPERATORS
	################################################################################
	
	updateBetaInitScale = 1
	# Scale operator for proportional shrinking-expanding move to update
	# initial phenotypic rate for rate regimes
	
	updateBetaShiftScale = 1
	# Scale operator for sliding window move to update initial phenotypic rate
	
	updateNodeStateScale = 1
	# Scale operator for sliding window move to update ancestral states
	# at internal nodes
	
	updateEventLocationScale = 0.05
	# Scale parameter for updating LOCAL moves of events on the tree
	# This defines the width of the sliding window proposal
	 
	updateEventRateScale = 4.0
	# Scale parameter (proportional shrinking/expanding) for updating
	# the rate parameter of the Poisson process 
	
	
	################################################################################
	# OPERATORS: MCMC MOVE FREQUENCIES
	################################################################################
	
	updateRateEventNumber = 1
	# Relative frequency of MCMC moves that change the number of events
	
	updateRateEventPosition = 1
	# Relative frequency of MCMC moves that change the location of an event
	# on the tree
	
	updateRateEventRate = 1
	# Relative frequency of MCMC moves that change the rate at which events occur 
	
	updateRateBeta0 = 1
	# Relative frequency of MCMC moves that change the initial phenotypic rate
	# associated with an event
	
	updateRateBetaShift = 1
	# Relative frequency of MCMC moves that change the exponential shift parameter
	# of the phenotypic rate associated with an event
	
	updateRateNodeState = 25
	# Relative frequency of MCMC moves that update the value of ancestral
	# character states. You have as many ancestral states as you have
	# internal nodes in your tree, so there are a lot of parameters:
	# you should update this much more often than you update the event-associated
	# parameters.
	
	updateRateBetaTimeMode = 0
	# Relative frequency of MCMC moves that flip the time mode
	# (time-constant <=> time-variable)
	
	localGlobalMoveRatio = 10.0
	# Ratio of local to global moves of events 
	
	
	################################################################################
	# INITIAL PARAMETER VALUES
	################################################################################
	
	betaInit = 0.5
	# Initial value of the phenotypic evolutionary process at the root of the tree
	
	betaShiftInit = 0
	# Initial value of the exponential change parameter for the phenotypic
	# evolutionary process at the root of the tree. A value of zero implies
	# time-constant rates
	
	initialNumberEvents = 0
	# Initial number of non-root processes
	
	
	################################################################################
	# METROPOLIS COUPLED MCMC
	################################################################################
	
	numberOfChains = 4
	# Number of Markov chains to run
	
	deltaT = 0.1
	# Temperature increment parameter. This value should be > 0
	# The temperature for the i-th chain is calculated as 1 / [1 + deltaT * (i - 1)]
	
	swapPeriod = 1000
	# Number of generations in which to propose a chain swap
	
	chainSwapFileName = chain_swap.txt
	# File name in which to output data about each chain swap proposal.
	# The format of each line is [generation],[rank_1],[rank_2],[swap_accepted]
	# where [generation] is the generation in which the swap proposal was made,
	# [rank_1] and [rank_2] are the chains that were chosen, and [swap_accepted] is
	# whether the swap was made. The cold chain has a rank of 1.')
	
	templates <- lapply(templates, function(x) strsplit(x, "\n")[[1]]);
	templates <- lapply(templates, function(x) gsub("\t", "", x));
	
	#identify appropriate template
	if (type == "diversification") {
		template <- templates$diversification;
	} else if (type == "trait") {
		template <- templates$trait;
	} else {
		stop("type must be either diversification or trait.");
	}
		
	#replace defaults with user-specified parameter values
	if (!is.null(params)) {
		params <- lapply(params, as.character);
		for (i in 1:length(params)) {
			paramName <- names(params)[i];
			paramName <- paste(paramName, " = ", sep='');
			if (!any(grepl(paramName, template))) {
				stop(paste(names(params)[i], " parameter not found in template.", sep = ""));
			} else {
				ind <- which(sapply(template, function(x) grepl(paramName, x)) == TRUE);
				
				#if multiple lines contain file name, 1 is likely correct, and the others are likely in comments
				if (length(ind) > 1) {
					isComment <- sapply(ind, function(x) grepl("#", template[x]), USE.NAMES = FALSE);
					ind <- setdiff(ind, ind[which(isComment == TRUE)]);
				}
				
				#if multiple matches, this should be fixed, so report
				if (length(ind) > 1 | length(ind) == 0) {
					stop(paste(names(params)[i], " returning multiple hits.", sep=""));
				}
				
				template[ind] <- gsub("=.+$", paste("= ", params[i], sep=""), template[ind]);
				
				#if parameter is commented out, uncomment it
				template[ind] <- gsub("# ", "", template[ind]);
			}
		}
	}
		
	#write file to disk
	write(template, file = file);
		
}
