#' @title 
#' Disease Model without Storage
#' 
#' @description 
#' This function is used to model disease (freedom and prevalence estimation) using 
#'     a Bayesian hierarchical model with one, two, or three levels of sampling. 
#'     The most general, the three level, would be represented by the hierarchy: 
#'     (region -> subzone -> cluster -> subject) and the two- and one- level sampling would 
#'     simply be the same with either subzone removed or subzone \emph{and} cluster 
#'     removed, respectively. This function is the 'no storage' model and hence stores 
#'     only immediately relevant output.
#'
#' @param H Number of subzones/states. Should be a 1 if there are only one or two levels 
#'     of sampling. Integer scalar.
#' @param k Number of clusters / farms / ponds / herds. Should be a 1 if there is only one
#'     level of sampling. Integer vector (\code{H} x 1).
#' @param n Number of subjects / animals / mussels / pigs per cluster (can differ among
#'     clusters). Integer vector (\code{sum(k)} x 1).
#' @param seasons Numeric season for  each cluster in the order: Summer (1),  Fall (2),
#'     Winter (3),  Spring (4). Integer vector (\code{sum(k)} x 1).
#' @param reps Number of (simulated) replicated data sets. Integer scalar.
#' @param MCMCreps Number of iterations in the MCMC chain per replicated data set. 
#'     Integer scalar.
#' @param poi The p(arameter) o(f) i(nterest) specifies one of the subzone-level prevalence
#'     (\code{gam}) or the cluster-level prevalence (\code{tau}), indicating which variable
#'     with which to compute the simulation output \code{p2.tilde}, \code{p4.tilde}, and
#'     \code{p6.tilde}. Character scalar.
#' @param y An optional input of sums of positive diagnostic testing results if one has a 
#'     specific set of diagnostic testing outcomes for  every subject (will simulate these 
#'     if  this is left as \code{NULL}). Integer matrix (\code{reps} x \code{sum(k)}).
#' @param mumodes Modes and (a) 95th percentiles for  mode \eqn{\leq}  0.50 or (b) 5th 
#'     percentiles for  mode \eqn{>} 0.5 for  season-specific mean prevalences for  
#'     diseased clusters in the order: Summer,  Fall,  Winter,  Spring. 
#'     Real matrix (4 x 2).
#' @param pi.thresh Threshold that we must show subject-level prevalence is below to 
#'     declare disease freedom. Real scalar.
#' @param tau.thresh Threshold that we must show cluster-level prevalence is below to  
#'     declare disease freedom. Real scalar.
#' @param gam.thresh Threshold that we must show subzone-level prevalence is below to  
#'     declare disease freedom. Real scalar.
#' @param tau.T Assumed true cluster-level prevalence (used to simulate data to feed 
#'     into the Bayesian model). Real scalar.
#' @param poi.lb,poi.ub Lower and upper bounds for posterior \code{poi} prevalences 
#'     to show ability to capture \code{poi} with certain probability. Real scalars.
#' @param p1 Probability we must show \code{poi} prevalence is below / above the thresholds 
#'     \code{pi.thresh} or \code{tau.thresh} or within specified bounds. Real scalar.
#' @param psi (Inversely related to) the variability of the subject-level prevalences 
#'     in diseased clusters. Real scalar. 
#' @param omegaparm Prior parameters for the beta-distributed input \code{omega}, which 
#'     is the probability that the disease is in the region. Real vector (2 x 1).
#' @param gamparm Prior parameters for the beta-distributed input \code{gam}, which 
#'     is the subzone-level prevalence (or the prevalence among subzones). Real 
#'     vector (2 x 1).
#' @param tauparm Prior parameters for the beta-distributed input \code{tau}, which is the
#'      cluster-level prevalence (or the prevalence among clusters). Real vector (2 x 1).
#' @param etaparm Prior parameters for the beta-distributed input \code{eta}, which is the 
#'     sensitivity of the diagnostic test. Real vector (2 x 1).
#' @param thetaparm Prior parameters for the beta-distributed input \code{theta}, which is 
#'     the specificity of the diagnostic test. Real vector (2 x 1).
#' @param burnin Number of MCMC iterations to discard from the beginning of the chain. 
#'     Integer scalar.
#'
#' @details
#' Note that this function performs in the same manner as \code{\link{EpiBayes_s}} except
#'     it doesn't store the posterior distributions of all of the parameters in the model
#'     and hence is a bit quicker to run.
#'     This model is a Bayesian hierarchical model that serves two main purposes: 
#'     \itemize{
#'         \item Simulation model: can simulate data under user-specified conditions 
#'         and run replicated data sets under the Bayesian model to observed the behavior 
#'         of the system under random realizations of simulated data. 
#'         \item Posterior inference model: can use actual observed data from the field, 
#'         run it through the Bayesian model,  and make inference on parameter(s) of
#'         interest using the posterior distribution(s).
#'     }
#'
#'     The posterior distributions are avaialable for a particular parameter, say tau, by
#'     typing \code{name_of_your_model$taumat}. 
#'     Note: be careful about the size of the \code{taumat} matrix you are calling. The 
#'     last index of any of the variables from above is the MCMC replications and so we 
#'     would typically always omit the last index when looking at any particular variable.
#'     \itemize{
#'         \item If  we want to look at the posterior distribution of the cluster level
#'         prevalence (\code{tau}) for  the first replication in the first subzone,  we 
#'         will note that taumat is a matrix with the first dimension indexed by 
#'         replication, the second dimension indexed by subzone, and the third dimension
#'         indexed by MCMC replications. Then, we will type something like 
#'         \code{name_of_your_model$taumat[1, 1, ]} to visually
#'         inspect the posterior distribution in the form of a vector. For the second
#'         replication in the fourth subzone, we can type 
#'         \code{name_of_your_model$taumat[2, 4, ]},  and so forth.
#'         Then,  we can make histograms of these distributions if we so desire by the 
#'         following code: 
#'         \code{hist(name_of_your_model$taumat[1, 1, ], col = "cyan");box("plot")}. 
#'         To observe a trace plot,  we can type: 
#'         \code{plot(name_of_your_model$taumat[1, 1, ], type = "l")} 
#'         for  all of the MCMC replications and we can look 
#'         at the trace plot after a burnin of 1000 iterations by typing: 
#'         \code{plot(name_of_your_model$taumat[1, 1, -c(1:1000)], type ="l")}.
#'     }
#'
#' @section Parameter Meanings:
#' Below we describe the meanings of the parameters in epidemiological language for ease 
#'     of implementation and elicitation of inputs to this function. 
#'
#' \tabular{lll}{
#'     Parameter \tab (3-level) Description \tab (2-level) Description \cr
#' 		 \code{omega} \tab Probability of disease being in the region. \tab Not used. \cr
#'       \code{gam} \tab Subzone-level (between-subzone) prevalence. \tab Probability of disease being in the region. \cr
#'       \code{z.gam} \tab Subzone-level (between-subzone) prevalence latent indicator variable. \tab Not used. \cr
#'		 \code{tau} \tab Cluster-level (between-cluster) prevalence. \tab Same as (3-level). \cr
#'		 \code{z.tau} \tab Cluster-level (between-cluster) prevalence latent indicator variable. \tab Same as (3-level). \cr
#'		 \code{pi} \tab Subject-level (within-cluster) prevalence. \tab Same as (3-level). \cr
#'		 \code{z.pi} \tab Subject-level (within-cluster) prevalence latent indicator variable. \tab Same as (3-level). \cr
#'		 \code{mu} \tab Mean prevalence among infected clusters. \tab Same as (3-level). \cr
#'		 \code{psi} \tab (Related to) variability of prevalence among infected clusters (inversely related so higher \code{psi} -> lower variance of prevalences among diseased clusters). \tab Same as (3-level). \cr
#'		 \code{eta} \tab Diagnostic test sensitivity. \tab Same as (3-level). \cr
#'		 \code{theta} \tab Diagnostic test specificity. \tab Same as (3-level). \cr
#'		 \code{c1} \tab Latent count of true positive diagnostic test results. \tab Same as (3-level). \cr
#'		 \code{c2} \tab Latent count of true negative diagnostic test results. \tab Same as (3-level). \cr
#' }
#'
#'     Note that in the code, each of these variables are denoted, for example, 
#'     \code{taumat} instead of plain \code{tau} to denote that that variable is a 
#'     matrix (or, more generally an array) of values.
#' 
#' @return
#' This function is the 'no storage' model meaning it does not return as many values as the
#'     'storage' model, \code{\link{EpiBayes_s}}. It returns enough to perform inference
#'     but \emph{not} to perform diagnostics for the model fit nor MCMC convergence. 
#'     If this is the first time running the model, we recommend that the user utilizes 
#'     the function \code{\link{EpiBayes_ns}} and diagnose issues before continuing with 
#'     this 'no storage' model. 
#'     Nevertheless, below are the outputs of this model.
#'
#' \tabular{lll}{
#'     \code{p2.tilde} \tab Real scalar \tab Proportion of simulated data sets that result in the probability of \code{poi} prevalence \emph{below} \code{poi.thresh} with probability \code{p1} \cr
#'     \code{p4.tilde} \tab Real scalar \tab Proportion of simulated data sets that result in the probability of \code{poi} prevalence \emph{above} \code{poi.thresh} with probability \code{p1} \cr
#'     \code{p6.tilde} \tab Real scalar \tab Proportion of simulated data sets that result in the probability of \code{poi} prevalence \emph{between} \code{poi.lb} and \code{poi.ub} with probability \code{p1} \cr
#'     \code{taumat} \tab Real array (\code{reps} x \code{H} x \code{MCMCreps}) \tab Posterior distributions of the cluster-level prevalence for all simulated data sets (i.e., \code{reps}) \cr
#'     \code{gammat} \tab  Real matrix (\code{reps} x \code{MCMCreps}) \tab Posterior distribution of the subzone-level prevalence \cr
#'     \code{ForOthers} \tab \tab Various other data not intended to be used by the user, but used to pass information on to the \code{plot}, \code{summary}, and \code{print} methods \cr
#' } 
#' 
#'     The posterior distribution output, say \code{taumat}, can be manipulated after it is
#'     returned with the \pkg{coda} package after it is converted to an \code{mcmc} object. 
#'
#' @references
#' Branscum, A., Johnson, W., and Gardner, I. (2006) Sample size calculations for disease 
#'     freedom and prevalence estimation surveys. \emph{Statistics in Medicine 25}, 
#'     2658-2674.
#'
#' @seealso
#' The function \code{\link{EpiBayes_s}} stores more output so the user may check 
#'     diagnostics, but is understandably slower to execute.
#'
#' @examples
#' testrun_nostorage = EpiBayes_ns(
#'		H = 2, 
#'		k = rep(30, 2),  
#'		n = rep(rep(150, 30), 2), 
#'		seasons = rep(c(1, 2, 3, 4), each = 15), 
#'		reps = 10, 
#'		MCMCreps = 10,
#'		poi = "tau",
#'		y = NULL,
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
#'		omegaparm = c(100, 1),
#'		gamparm = c(100, 1), 
#'		tauparm = c(1, 1), 
#'		etaparm = c(100, 6), 
#'		thetaparm = c(100, 6),
#'		burnin = 1 
#'		)
#'
#' testrun_nostorage
#' print(testrun_nostorage)
#' testrun_nostoragesummary = summary(testrun_nostorage, prob = 0.90, n.output = 5)
#' testrun_nostoragesummary
#' plot(testrun_nostorage$taumat[1, 1, ], type = "l")
#' plot(testrun_nostorage$gammat[1, ], type = "l")
#' 
#' ## Can't look at any other posterior distributions other than gam and tau
#' \dontrun{
#'     plot(testrun_nostorage$pimat[1, 1, ], type = "l")
#' 	   plot(testrun_nostorage$omegamat[1, ], type = "l")
#' }
#' 
#' @export
#' @name EpiBayes_ns
#' @import compiler
#' @import epiR
EpiBayes_ns = function(
					  H, 
					  k, 
					  n, 
					  seasons, 
					  reps, 
					  MCMCreps,
					  poi = "tau",
					  y = NULL,
					  mumodes = matrix(c(
							0.50, 0.70, 
							0.50, 0.70, 
							0.02, 0.50, 
							0.02, 0.50
							), 4, 2, byrow = TRUE
						),
					  pi.thresh = 0.02,
					  tau.thresh = 0.05, 
					  gam.thresh = 0.10,
					  tau.T = 0, 
					  poi.lb = 0, 
					  poi.ub = 1, 
					  p1 = 0.95, 
					  psi = 4,
					  omegaparm = c(100, 1),
					  gamparm = c(100, 1), 
					  tauparm = c(1, 1), 
					  etaparm = c(100, 6), 
					  thetaparm = c(100, 6),
					  burnin = 1000 
					  ){	
	##Error Handling
		#Error handling making sure we have all subzone sizes
		if (H !=  length(k))	stop("DIMENSION MISMATCH: The number of subzones in input 'H' is not equal to the number of subzones by the length of input vector 'k'.")
	
		#Error handling making sure we have all subzone sizes
		if (length(n) !=  sum(k))	stop("DIMENSION MISMATCH: The number of clusters provided by the length of input vector 'n' is not equal to the number of clusters provided by input 'sum(k)'.")
		#Error handling making sure we have mean prevalences among diseased clusters for  all four seasons in order of Summer,  Fall,  Winter,  Spring
		if (nrow(mumodes) !=  4) stop("DIMENSION MISMATCH: The number of mean prevalences and their 95th percentiles ('mumodes') is not equal to 4. Make sure to provide one for  each season in order of: Summer,  Fall,  Winter,  Spring")
		#Similar,  checking number of seasons given
		if (length(seasons) !=  sum(k)) stop("DIMENSION MISMATCH: The 'seasons' input should be a vector of length 'sum(k)'. Give every cluster its own numeric season as: Summer (1),  Fall (2),  Winter (3),  Spring (4)")
		# Checking if y is a matrix
		if (!is.null(y) & !is.matrix(y)) stop("OBJECT TYPE ERROR: The input matrix 'y' should be a matrix with 'reps' rows and 'sum(k)' columns. Or, it should be 'NULL'.")
		# Checking size of matrix
		if (is.matrix(y)){
			if(!isTRUE(all.equal(dim(y), c(reps, sum(k))))) stop("DIMENSION MISMATCH: The input matrix 'y' should have 'reps' rows and 'sum(k)' columns. Or, it should be 'NULL'.")	
		}
		
	## Declare default point mass parameters that define a beta distribution that looks like a point mass at zero
	pmparm = c(1, 100)
	
	## Declare number of seasons
	n.seasons = nrow(mumodes)	
	
	## Store input y matrix
	iny = y
	
	## Create empty vectors for  simulation perfor mance outputs
	p2.tilde = numeric(reps)
	p4.tilde = numeric(reps)
	p6.tilde = numeric(reps)
	
	## Create empty array for tau and gamma (tracking this parameter only, so outside of the reps loop)
	taumat = array(0, c(reps, H, MCMCreps))
	gammat = matrix(0, reps, MCMCreps)
	
	##Initial values for tau and gamma
	taumat[, , 1] = rbeta(reps * H, tauparm[1], tauparm[2])
	gammat[, 1] = rbeta(reps, gamparm[1], gamparm[2])

	##Outer loop for  replicated data sets
	for (r in 1:reps){
		
		#Set up empty arrays to hold realizations
		z.pimat = matrix(0, sum(k), MCMCreps)
			z.pimat2 = matrix(0, sum(k), MCMCreps)
		pimat = matrix(0, sum(k), MCMCreps)
		c1mat = matrix(0, sum(k), MCMCreps)
		c2mat = matrix(0, sum(k), MCMCreps)
		z.taumat = matrix(0, H, MCMCreps)
			z.taumat2 = matrix(0, H, MCMCreps)
		z.gammat = numeric(MCMCreps)
			z.gammat2 = numeric(MCMCreps)
		omegamat = numeric(MCMCreps)
		etamat = numeric(MCMCreps)
		thetamat = numeric(MCMCreps)
		mumat = matrix(0, n.seasons, MCMCreps)
		
		#Set up and populate matrix of prior parameters for  season - specif ic mu values
		muparm = matrix(0, nrow = n.seasons, ncol = 2)
		midflag = (mumodes[, 1]<= 0.5)
		for (i in 1:nrow(mumodes)){
				mu.epi = epiR::epi.betabuster(mumodes[i, 1], 0.95, 1 - midflag[i], mumodes[i, 2])
				muparm[i, ] = c(mu.epi$shape1, mu.epi$shape2)
			}
			
		#Generate some sample data under prior distributions and tau.T
		if (!is.null(iny)){
			y = iny[r, ]
			} else {
			#Null sensitivitites and specif icities under tau.T
			eta.0 = rbeta(sum(k), etaparm[1], etaparm[2])
			theta.0 = rbeta(sum(k), thetaparm[1], thetaparm[2])
			
			#Null subject - level prevalences under tau.T		
			mu.0 = rbeta(sum(k), rep(muparm[, 1], table(factor(seasons, levels = c(1, 2, 3, 4)))), rep(muparm[, 2], table(factor(seasons, levels = c(1, 2, 3, 4)))))
			psi.0 = psi
			z.pi.0 = rbinom(sum(k), 1, tau.T)
			pi.0 = rbeta(sum(k), mu.0 * psi.0, (1 - mu.0) * psi.0) * z.pi.0
			
			#Simulate the outcomes under prior distributions and tau.T
			y = rbinom(sum(k), n, pi.0 * eta.0 + (1 - pi.0) * (1 - theta.0))
			} 
		
		#Set initial values from prior distributions under tau.T
		omegamat[1] = rbeta(1, omegaparm[1], omegaparm[2])
		z.gammat[1] = rbinom(1, 1, omegamat[1])
		gammat[r, 1] = gammat[r, 1] * z.gammat[1]
			z.gammat2[1] = ifelse(gammat[r, 1] >= gam.thresh, 1, 0)
		z.taumat[, 1] = rbinom(H, 1, gammat[r, 1])
		taumat[r, , 1] = taumat[r, , 1] * z.taumat[, 1]
			z.taumat2[, 1] = ifelse(taumat[r, , 1] >= tau.thresh, 1, 0)
		mumat[, 1] = rbeta(n.seasons, rep(muparm[, 1], table(factor(seasons, levels = c(1, 2, 3, 4)))), rep(muparm[, 2], table(factor(seasons, levels = c(1, 2, 3, 4)))))
		z.pimat[, 1] = rbinom(sum(k), 1, taumat[, , 1]) 
		pimat[, 1] = rbeta(sum(k), mumat[, 1] * psi, (1 - mumat[, 1]) * psi) * z.pimat[, 1]
			z.pimat2[, 1] = ifelse(pimat[, 1] >= pi.thresh, 1, 0)
		etamat[1] = rbeta(1, etaparm[1], etaparm[2])
		thetamat[1] = rbeta(1, thetaparm[1], thetaparm[2])
		c1mat[, 1] = rbinom(sum(k), y, thetamat[1])
		c2mat[, 1] = rbinom(sum(k), n - y, etamat[1])

		#Sampling from conditionals outer loop
		for (t in 2:MCMCreps){
		
			#Sampling from conditionals for  parameters which have unique values per subzone
			for (h in 1:H){
				
				#Sampling from conditionals for  parmameters which have unique values per cluster
				for (i in 1:k[h]){
					#Latent pi indicator
						#Get the components of the probability of declaring a cluster diseased
						winclust.ind = sum(k[0:(h - 1)]) + i
						
						piprobbit.a = (taumat[r, h, t - 1] * dbeta(pimat[winclust.ind, t - 1], mumat[seasons[winclust.ind], t - 1] * psi, (1 - mumat[seasons[winclust.ind], t - 1]) * psi))
						piprobbit.b = ((1 - taumat[r, h, t - 1]) * dbeta(pimat[winclust.ind, t - 1], pmparm[1], pmparm[2]))
							#To control for  pi = 0 case (assume the exactly zero values come from the point mass)
							if (is.infinite(piprobbit.a) | is.nan(piprobbit.a)){piprobbit.a = 0} 
							#To control for  both zero case (assume these x^ - 98 values come from the point mass)
							if (piprobbit.a == 0 & piprobbit.b == 0){piprobbit.b = 1} 
					z.pimat[winclust.ind, t] = rmultinom(1, 1, c(piprobbit.a, piprobbit.b))[1]  # ifelse(pimat[winclust.ind, t - 1] >= pi.thresh, 1, 0) # rmultinom(1, 1, c(piprobbit.a, piprobbit.b))[1] 
					
					#Pi
					pimat[winclust.ind, t] = rbeta(1, z.pimat[winclust.ind, t] * (mumat[seasons[winclust.ind], t - 1] * psi - 1) + (1 - z.pimat[winclust.ind, t]) * (pmparm[1] - 1) + c1mat[winclust.ind, t - 1] + n[winclust.ind] - y[winclust.ind] - c2mat[winclust.ind, t - 1] + 1,  z.pimat[winclust.ind, t] * ((1 - mumat[seasons[winclust.ind], t - 1]) * psi - 1) + (1 - z.pimat[winclust.ind, t]) * (pmparm[2] - 1) + y[winclust.ind] - c1mat[winclust.ind, t - 1] + c2mat[winclust.ind, t - 1] + 1)
						#In order to avoid the problem of setting mu target to 0 if  a pi value is equal to 1 exactly
						if (pimat[winclust.ind, t] == 1){pimat[winclust.ind, t] = 0.999999} 
						#Forces "essentially" zero prevalence estimates to exactly zero to better determine which prevalences come from which mixture component while still allowing for  jumps between the mixture components if  the sample data are pursuasive enough
						if (round(pimat[winclust.ind, t], 80) == 0){pimat[winclust.ind, t] = 0} 
					# Calculate latent variable in terms of design prevalence to pass to higher levels	
					z.pimat2[winclust.ind, t] = ifelse(pimat[winclust.ind, t - 1] >= pi.thresh, 1, 0)
												
					#True positives
					c1mat[winclust.ind, t] = rbinom(1, y[winclust.ind],  (pimat[winclust.ind, t] * etamat[t - 1]) / (pimat[winclust.ind, t] * etamat[t - 1] + (1 - pimat[winclust.ind, t]) * (1 - thetamat[t - 1])))
					
					#True negatives
					c2mat[winclust.ind, t] = rbinom(1, n[winclust.ind] - y[winclust.ind],  ((1 - pimat[winclust.ind, t]) * thetamat[t - 1]) / ((1 - pimat[winclust.ind, t]) * thetamat[t - 1] + pimat[winclust.ind, t] * (1 - etamat[t - 1])))
					
					}
					
			#Latent tau indicator
				#Get the components of the probability of declaring a region diseased
				tauprobbit.a = (gammat[r, t - 1] * dbeta(taumat[r, h, t - 1], tauparm[1], tauparm[2]))
				tauprobbit.b = ((1 - gammat[r, t - 1]) * dbeta(taumat[r, h, t - 1], pmparm[1], pmparm[2]))
					#To control for  pi = 0 case (assume the exactly zero values come from the point mass)
					if (is.infinite(tauprobbit.a) | is.nan(tauprobbit.a)){tauprobbit.a = 0} 
					#To control for  both zero case (assume these x^ - 98 values come from the point mass)
					if (tauprobbit.a == 0 & tauprobbit.b == 0){tauprobbit.b = 1} 
			z.taumat[h, t] = rmultinom(1, 1, c(tauprobbit.a, tauprobbit.b))[1] # ifelse(taumat[r, h, t - 1] >= tau.thresh, 1, 0) # rmultinom(1, 1, c(tauprobbit.a, tauprobbit.b))[1] 
			 
			#Tau
			taumat[r, h, t] = rbeta(1, sum(z.pimat2[(winclust.ind - i + 1):winclust.ind, t]) + z.taumat[h, t] * (tauparm[1] - 1) + (1 - z.taumat[h, t]) * (pmparm[1] - 1) + 1,  k[h] - sum(z.pimat2[(winclust.ind - i + 1):winclust.ind, t]) + z.taumat[h, t] * (tauparm[2] - 1) + (1 - z.taumat[h, t]) * (pmparm[2] - 1) + 1)
					#Forces "essentially" zero prevalence estimates to exactly zero to better determine which prevalences come from which mixture component while still allowing for  jumps between the mixture components if  the sample data are pursuasive enough
					if (round(taumat[r, h, t], 80) == 0){taumat[r, h, t] = 0}
					
			# Calculate latent variable in terms of design prevalence to pass to higher levels	
			z.taumat2[h, t] = ifelse(taumat[r, h, t - 1] >= tau.thresh, 1, 0)
				}
			
			##MH step for  mu			
				#Loop through the season - specif ic prevalences			
				for (s in 1:n.seasons){
		
					#Count the number of clusters declared diseased in each unique season
					whichzone = which(z.pimat2[, t] == 1 & seasons == s)
					n.whichzone = sum(z.pimat2[whichzone, t])
					
					#Set the current value of mu
					mumh.t = mumat[s, t - 1]
					
					#Estimate the mean and variance of the posterior using MOM estimation to simulate from proposal distribution
					mumh.hatt = sum(pimat[whichzone, t]) / n.whichzone
					mumh.varhatt = sum((pimat[whichzone, t]  -  mumh.hatt)^2) / n.whichzone^2
					if (n.whichzone>1){
						mumh.alphastar = ((mumh.hatt * (1 - mumh.hatt) / mumh.varhatt) - 1) * mumh.hatt
						mumh.betastar = ((mumh.hatt * (1 - mumh.hatt) / mumh.varhatt) - 1) * (1 - mumh.hatt)
					} else{
						mumh.alphastar = 0
						mumh.betastar = 0
						}
						
					#Simulate from proposal distribution
					mumh.star = rbeta(1, mumh.alphastar + muparm[s, 1] - 1, mumh.betastar + muparm[s, 2] - 1)
						# To avoid numerical issues when mumh.star = 1
						if(mumh.star == 1){mumh.star = 0.9999}
					
					#Decide whether or not to accept or reject proposed value
					mumh.ratio = (EpiBayes::utils_mumhtargetdist(mumh.star, psi, pimat[, t], z.pimat[, t], muparm[s, ]) * EpiBayes::utils_mumhproposaldist(mumh.t, muparm[s, 1], muparm[s, 2], mumh.alphastar, mumh.betastar)) / (EpiBayes::utils_mumhtargetdist(mumh.t, psi, pimat[, t], z.pimat[, t], muparm[s, ]) * EpiBayes::utils_mumhproposaldist(mumh.star, muparm[s, 1], muparm[s, 2], mumh.alphastar, mumh.betastar))
						#To protect against possibility of having both numerator and denominator of MH ratio zero
						if (is.nan(mumh.ratio)){
							mumh.ratio = 0
						}
					mumh.acceptp = min(mumh.ratio, 1)
					mumh.acceptbern = rbinom(1, 1, mumh.acceptp)
						if (mumh.acceptbern == 1){
							mumat[s, t] = mumh.star
						} else{
							mumat[s, t] = mumh.t
							}
					}
			
			#Latent gamma indicator
				#Get the components of the probability of declaring a region diseased
				gamprobbit.a = (omegamat[t - 1] * dbeta(gammat[r, t - 1], gamparm[1], gamparm[2]))
				gamprobbit.b = ((1 - omegamat[t - 1]) * dbeta(gammat[r, t - 1], pmparm[1], pmparm[2]))
					#To control for  pi = 0 case (assume the exactly zero values come from the point mass)
					if (is.infinite(gamprobbit.a) | is.nan(gamprobbit.a)){gamprobbit.a = 0}
					#To control for  both zero case (assume these x^ - 98 values come from the point mass)
					if (gamprobbit.a == 0 & gamprobbit.b == 0){gamprobbit.b = 1}  
			z.gammat[t] = rmultinom(1, 1, c(gamprobbit.a, gamprobbit.b))[1]
					
			#Gamma
			gammat[r, t] = rbeta(1, sum(z.taumat2[, t]) + z.gammat[t] * (gamparm[1] - 1) + (1 - z.gammat[t]) * (pmparm[1] - 1) + 1, H - sum(z.taumat2[, t]) + z.gammat[t] * (gamparm[2] - 1) + (1 - z.gammat[t]) * (pmparm[2] - 1) + 1)
			
			# Calculate latent variable in terms of design prevalence to pass to higher levels
			z.gammat2[t] = ifelse(gammat[r, t] >= gam.thresh, 1, 0)
			
			#Omega
			omegamat[t] = rbeta(1, z.gammat2[t] + omegaparm[1], 1 - z.gammat2[t] + omegaparm[2])
				
			#Eta
			etamat[t] = rbeta(1, sum(c1mat[, t]) + etaparm[1],  sum(n - y[] - c2mat[, t]) + etaparm[2])
			
			#Theta
			thetamat[t] = rbeta(1, sum(c2mat[, t]) + thetaparm[1],  sum(y - c1mat[, t]) + thetaparm[2])
												
		}
	
	##Calculate items to be saved for  post - analysis
	if (poi == "tau"){
			
		#Calculate p2,  the proportion of simulations with probability of disease below threshold with p1 probability (probability of disease - freedom)
		p2.tilde[r] = mean(taumat[r, , -c(1:burnin)] <=  tau.thresh)
	
		#Calculate p4,  the proportion of simulations with probability of disease above threshold with p1 probability (probability of disease)
		p4.tilde[r] = mean(taumat[r, , -c(1:burnin)] > tau.thresh)
	
		#Calculate cover2,  the proportion of simulations with probability mass within the specif ied interval for  tau at or above the specif ied level,  p1
		p6.tilde[r] = mean(taumat[r, , -c(1:burnin)] >=  poi.lb & taumat[r, , -c(1:burnin)] <=  poi.ub)
			
	}else if (poi == "gamma"){
		
		#Calculate p2,  the proportion of simulations with probability of disease below threshold with p1 probability (probability of disease - freedom)
		p2.tilde[r] = mean(gammat[r, -c(1:burnin)] <=  gam.thresh)
	
		#Calculate p4,  the proportion of simulations with probability of disease above threshold with p1 probability (probability of disease)
		p4.tilde[r] = mean(gammat[r, -c(1:burnin)] > gam.thresh)
	
		#Calculate cover2,  the proportion of simulations with probability mass within the specif ied interval for  tau at or above the specif ied level,  p1
		p6.tilde[r] = mean(gammat[r, -c(1:burnin)] >=  poi.lb & gammat[r, -c(1:burnin)] <=  poi.ub)
				
	} else{
		p2.tilde = NULL
		p4.tilde = NULL
		p6.tilde = NULL
			}
			
	}
	
	p2.tilde = mean(p2.tilde >=  p1)
	p4.tilde = mean(p4.tilde >=  p1)
	p6.tilde = mean(p6.tilde >=  p1)
	
	EpiBayes3ns.out = list(
						"p2.tilde" = p2.tilde, 
						"p4.tilde" = p4.tilde, 
						"p6.tilde" = p6.tilde, 
						"taumat" = taumat,
						"gammat" = gammat,
						"omegamat" = NULL,
						"z.gammat" = NULL,
						"z.taumat" = NULL,
						"pimat" = NULL, 
						"z.pimat" = NULL, 
						"mumat" = NULL,
						"psi" = NULL,
						"etamat" = NULL, 
						"thetamat" = NULL,
						"c1mat" = NULL, 
						"c2mat" = NULL, 
						"mumh.tracker" = NULL,
						"y" = NULL,
						"ForOthers" = list(
										"H" = H, 
										"k" = k, 
										"n" = n, 
										"seasons" = seasons, 
										"reps" = reps, 
										"MCMCreps" = MCMCreps, 
										"mumodes" = mumodes, 
										"tau.T" = tau.T, 
										"poi" = poi, 
										"pi.thresh" = pi.thresh,
										"tau.thresh" = tau.thresh,
										"gam.thresh" = gam.thresh, 
										"poi.lb" = poi.lb, 
										"poi.ub" = poi.ub, 
										"p1" = p1,
										"psi" = psi,
										"omegaparm" = omegaparm,
										"gamparm" = gamparm, 
										"tauparm" = tauparm, 
										"etaparm" = etaparm, 
										"thetaparm" = thetaparm, 
										"burnin" = burnin
										)
						)
										
	class(EpiBayes3ns.out) = "eb"									

	return(invisible(EpiBayes3ns.out))
	
}

#Compile the hard storage model
EpiBayes_ns = compiler::cmpfun(EpiBayes_ns)