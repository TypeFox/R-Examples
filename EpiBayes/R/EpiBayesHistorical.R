#' @title 
#' Historical Aggregation for Disease Model
#' 
#' @description
#' This function aggregates results concerning cluster-level prevalence from the
#'     disease model from function \code{\link{EpiBayes_ns}}
#'     across a given number of time periods. Generally, this is 
#'     accomplished by starting out with a user-defined prior distribution on the 
#'     cluster-level prevalence (\code{tau}), updating this prior using observed data,
#'     using this posterior distribution as the prior distribution in the second time
#'     period, and repeating this process for all time periods -- we hope to implement
#'     a way to incorporate a measurement of introduction risk of the disease in 
#'     between time periods, but for now we assume that the disease retains its original
#'     properties among all time periods.
#'
#' @param input.df Data frame of input values that \emph{must} be supplied by the user, and 
#'     will be passed to the function \code{\link{EpiBayes_ns}}. The matrix will have one 
#'     row per cluster and will have five columns. The columns should have the order: 
#'     \strong{Time Period}, \strong{Subzone}, \strong{Cluster Size}, \strong{Season}, and 
#'     \strong{Positive Diagnostic Test Results (y)}. See \code{Details}
#'     for more information. Real matrix (\code{sum(k)} x 5).
#' @param orig.tauparm The prior parameters for the beta-distributed cluster-level 
#'     prevalence we assume to hold before our first time period. Real vector (2 x 1).
#' @param MCMCreps Number of iterations in the MCMC chain per replicated data set. 
#'     Integer scalar.
#' @param burnin Number of MCMC iterations to discard from the beginning of the chain. 
#'     Integer scalar.
#' @param ... Additional arguments that will be passed to \code{\link{EpiBayes_ns}}. 
#'     Otherwise, the default values will be used.
#'
#' @return
#' The returned values are given in a list. They are as follows.
#'
#' \tabular{lll}{
#'     Output \tab Attributes \tab Description \cr
#'     \code{RawPost} \tab List: Length - (number of periods), Elements - Real arrays (\code{reps} x \code{H} x \code{MCMCreps}) \tab Posterior distributions for the cluster-level prevalences for each subzone from all time periods \cr
#'     \code{BetaBusterEst} \tab List: Length - (number of periods), Elements - Real vectors (2 x 1) \tab Estimated posterior distributions for the cluster-level prevalences for each subzone from all time periods using moment-matching to the closest beta distribution by the function \code{\link[epiR]{epi.betabuster}} \cr
#'     \code{ForOthers} \tab  \tab Various other data not intended to be used by the user, but used to pass information on to the \code{plot}, \code{summary}, and \code{print} methods \cr
#' }
#'
#' @details
#' The \code{input.df} should have the following columns, in this order:
#'     \itemize{
#'         \item \strong{Time Period}: vector of codes for unique time periods. For 
#'             example, could be a vector of periods: \code{c(2015, 2015, 2016, ...)}.
#'         \item \strong{Subzone}: vector of codes for unique subzones. Should be the same
#'             for all rows if using one- or two-level sampling. For example, could be a 
#'             vector of names of a particular subzone: \code{c("CO", "CO", "IN", ...)}.
#'         \item \strong{Cluster Size}: vector of integers denoting the number of subjects 
#'             within that particular cluster. For example: \code{c(100, 500, 250, ...)}.
#'         \item \strong{Season}: vector of codes for season in which observed data were 
#'             collected. Must adhere to the requrirement that (1) denotes Summer, 
#'             (2) Fall, (3) Winter, and (4) Spring. For example: \code{c(1, 1, 4, ...)}.
#'         \item \strong{Positive Diagnostic Test Results (y)}: vector of integers
#'             denoting the observed number of positive diagnostic test results
#'             within that particular cluster. For example: \code{c(0, 4, 1, ...)}. Note:
#'             if so desired, the user may let the model generate sample data automatically
#'             when there is no concrete sample data with which to work.
#'     }
#'
#' @examples
#' ## Construct input data frame with columns Year, Subzone, Cluster size, Season, and Number positives
#' year = rep(c("Period 1", "Period 2", "Period 3"), c(60, 60, 60))
#' subz = rep(rep(c("Subzone 1", "Subzone 2"), c(25, 35)), 3)
#' size = rep(100, 3 * 60)
#' season = rep(rep(c(1,2), each = 30), 3)
#' y = matrix(c(
#'     rep(10, 15), rep(0, 10),  # Period 1: Subzone 1
#'     rep(0, 35),  # Period 1: Subzone 2
#'     rep(10, 15), rep(0, 10),  # Period 2: Subzone 1
#'     rep(10, 10), rep(0, 25),  # Period 2: Subzone 2
#'     rep(25, 25), # Period 3: Subzone 1
#'     rep(25, 10), rep(0, 25)  # Period 3: Subzone 2
#'     ),
#'     ncol = 1
#' )
#' 
#' testrun_historical_inputdf = data.frame(year, subz, size, season, y)
#' 
#' testrun_historical = EpiBayesHistorical(
#'		input.df = testrun_historical_inputdf,
#'		orig.tauparm = c(1, 1),
#'		burnin = 1,
#'		MCMCreps = 5,
#'		poi = "tau",
#'		mumodes = matrix(c(
#'			0.50, 0.70, 
#'			0.50, 0.70, 
#'			0.02, 0.50, 
#'			0.02, 0.50
#'			), 4, 2, byrow = TRUE
#'		),
#'		pi.thresh = 0.05, 
#'	    tau.thresh = 0.02,
#'      gam.thresh = 0.10,
#'		tau.T = 0, 
#'		poi.lb = 0, 
#'		poi.ub = 1, 
#'		p1 = 0.95, 
#'		psi = 4,
#'		omegaparm = c(1, 1),
#'		gamparm = c(1, 1), 
#'		etaparm = c(10, 1), 
#'		thetaparm = c(10, 1)
#'		)
#'
#' testrun_historical
#' plot(testrun_historical)
#' testrun_historicalsummary = summary(testrun_historical, sumstat = "quantile", 
#'     prob = 0.99, time.labels = c("Period 1", "Period 2", "Period 3"))
#' testrun_historicalsummary
#' plot(testrun_historicalsummary)
#' 
#' @export
#' @name EpiBayesHistorical
#' @import compiler
#' @import epiR 
EpiBayesHistorical = function(
						input.df,
				 		orig.tauparm,
				 		# orig.gamparm = NULL,
				 		# intro.prob = NULL,
				  		MCMCreps,
				  		burnin = 1000,
			      		...
						){
						
	# Separate each individual vector of inputs into useful list to be fed into historical function from the 'nice' form user inputs data
	input.df = data.frame(input.df)
	names(input.df) = c("period", "subzone", "cluster.size", "season", "positive")
	input.df$period = as.numeric(input.df$period)
	input.df$subzone = factor(input.df$subzone)
 	input.df$cluster.size = as.integer(input.df$cluster.size)
 	input.df$season = as.numeric(input.df$season)
 	input.df$positive = as.numeric(input.df$positive)
	
	n.periods = length(unique(input.df$period))

	input.H = list()
	input.k = list()
	input.n = list()
	input.seasons = list()
	input.y = list()
	input.subnames = list()
	
	for(i in 1:n.periods){
		which.currentperiod = which(input.df$period == unique(input.df$period)[i])
		# H
		input.H[[i]] = length(unique(input.df$subzone[which.currentperiod]))
		# k	
			table.k = table(input.df$subzone[which.currentperiod])
		input.k[[i]] = table.k[which(table.k != 0)]
		# n
		input.n[[i]] = input.df$cluster.size[which.currentperiod]
		# seasons
		input.seasons[[i]] = input.df$season[which.currentperiod]
		# all zero y matrices
		input.y[[i]] = matrix(input.df$positive[which.currentperiod], 1, sum(input.k[[i]]))
		# subzones
		input.subnames[[i]] = names(input.k[[i]])
	}
	
	# Alphabetize unique subzone names if not already
	unique.subzones = sort(unique(input.df$subzone))
					
	# Store number of subzones
	n.subzones = length(unique.subzones)
							  	
	# Initialize some variable to keep track of the posterior distributions among periods
	periodi.run = list()
	periodi.tau = list()
	periodi.betafit = list()
	periodi.betabust = list()

	# Initialize vector for holding names for legend
	periodi.names = character(n.periods)
	
	# Initialize variable to store max density for plotting purposes
	maxdens = matrix(0, n.periods, n.subzones)

	# Initialize tau parameters with original parameters supplied by the user
	current.tauparm = orig.tauparm
	
	# # Initialize gam parameters with original parameters supplied by the user
	# if (!is.null(orig.gamparm)){
		# current.gamparm = orig.gamparm
		# useorig = TRUE
		# } else if (!is.null(list(...)$gamparm)){
			# current.gamparm = list(...)$gamparm
			# useorig = FALSE
		# } else{
			# current.gamparm = c(1000, 1)
			# useorig = FALSE
			# }
				
	# Begin for loop for periods
	for (i in 1:n.periods){ 
		
		# Run the Bayesian model on each periods' parmeters, the optional additional parameters, and the current tau parameters
		periodi.run[[i]] = EpiBayes::EpiBayes_ns(
									H = input.H[[i]],
									k = input.k[[i]],
									n = input.n[[i]],
									seasons = input.seasons[[i]],
									y = input.y[[i]],
									tauparm = current.tauparm,
									MCMCreps = MCMCreps,
									burnin = burnin,
									reps = 1,
									# gamparm = current.gamparm,
									...
									)

		# Store output tau chains for each subzone
			# Create empty matrix for tau chains
			periodi.tau[[i]] = matrix(0, n.subzones, length(periodi.run[[i]]$taumat[1, 1, ]))
			# Store tau chains in appropriate spaces
			periodi.tau[[i]][c(match(input.subnames[[i]], unique.subzones)), ] = periodi.run[[i]]$taumat[1, , ]

		# Create empty matrix for BetaBuster output
		periodi.betabust[[i]] = matrix(0, n.subzones, 2)
		
		for(j in match(input.subnames[[i]], unique.subzones)){
			# Store the density smooth estimate in a list
			periodi.betafit = density(periodi.tau[[i]][j ,-c(1:burnin)], from = 0, to = 1)
			
			# Store the mode, 95% percentile, and use BetaBuster to estimate the best-fitting beta distribution to the smoothed estimate of the latest posterior distribution of tau
				# periodi.mode = periodi.betafit$x[which(periodi.betafit$y == max(periodi.betafit$y))]
				# periodi.95per = quantile(periodi.tau[[i]][j, -c(1:burnin)], 0.95)
				temp.betabust = EpiBayes::utils_newalphbet(mean(periodi.tau[[i]][j, -c(1:burnin)]), var(periodi.tau[[i]][j, -c(1:burnin)]))
				
			# Store the BetaBuster fit estimated parameters in a matrix
			periodi.betabust[[i]][j, 1] = temp.betabust[1]
			periodi.betabust[[i]][j, 2] = temp.betabust[2]
			
			# Store the maximum value of the density during the current period
			maxdens[i, j] = max(periodi.betafit$y)
		}

		# Redefine current.tauparm to reflect update
		current.tauparm = apply(matrix(periodi.betabust[[i]][match(input.subnames[[i]], unique.subzones),], ncol = 2), 2, mean)
		
		# # Redefine current.gamparm to reflect update if prior intro is supplied
		# if (useorig){
			
		# }
				
		# Store the name of the period for legend plotting
		periodi.names[i] = paste("Period ", i)
		
	}
	
	EpiBayesHistorical.out = list(
									"RawPost" = periodi.tau, 
									"BetaBusterEst" = periodi.betabust, 
									"ForOthers" = list(
													n.periods = n.periods, 
											 		orig.tauparm = orig.tauparm,
											  		burnin = burnin,
											  		input.df = input.df, 
													unique.subzones = unique.subzones, 
													maxdens = maxdens 
													)
									)
	
	class(EpiBayesHistorical.out) = "ebhistorical"
	
	return(invisible(EpiBayesHistorical.out))
}