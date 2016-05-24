#' Multivariate Normal Order-statistics Model.
#'
#' Using MCMC methods to fit the MVNOS model. Please install JAGS 3.X (\url{http://mcmc-jags.sourceforge.net}) and rjags (\url{https://cran.r-project.org/package=rjags}) at first.
#'
#' @param y 	:an n*k matrix, observed data, each row is an individual's rank of items
#' @param p 	:number of parameters in MVNOS model
#' @param Z 	:a n*k*p array of covariates associated with all judges
#' @param beta0 :a 1*p matrix, prior normal distribution mean parameters
#' @param A0 	:a p*p matrix, prior normal distribution variance-covariance matrix
#' @param alpha :scalar, prior Wishart distribution degree of freedom
#' @param P 	:a (k-1)*(k-1) matrix, prior Wishart distribution scale matrix
#' @param BURN_IN_ITERATIONS	:number of iterations to burn-in at first 
#' @param MAX_ITERATIONS	:full sample iterations
#' @param DRAW_CYCLE	:reduce the full sample by draw-cycle(e.g. draw every 20th draw from the full sample)
#' @return A list of Gibbs sampling traces
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' # APA data application
#' # It will take about 10 minutes to run the demo.
#' data(APA)
#' y = freq2case(APA, freq.col = 1)
#' y = 6 - y
#' # number of observed judges
#' n = dim(y)[1]
#' # number of items
#' k = dim(y)[2] 
#' # number of parameteros of beta
#' p = k 
#' beta0 = rep(0, p)
#' alpha = k + 1 
#' A0 = diag(100, ncol = p, nrow = p)
#' P = diag(k + 1, ncol = k - 1, nrow = k - 1)
#' # Construct Z
#' Z = array(0, dim = c(n, k, p))
#' for (j in 1:n)
#' {
#'    Z[j, , ] = diag(1, nrow= k, ncol = p)
#' }
#' # Total iterations of Gibbs sampling
#' MAX_ITERATIONS = 10000
#' # Number of iterations to be reduced(burnt in)
#' BURN_IN_ITERATIONS = 1000
#' # Run the model, time consuming
#' # output_list = mvnos.model(y, p, Z, beta0, A0, alpha, P, 
#' # MAX_ITERATIONS = MAX_ITERATIONS, BURN_IN_ITERATIONS = BURN_IN_ITERATIONS)
#' @references Yu, P. L. H. (2000). Bayesian analysis of order-statistics models for ranking data. Psychometrika, 65(3):281-299.

mvnos.model <- function(y, p, Z, beta0 = NULL, A0 = NULL, alpha = NULL, P = NULL, 
	BURN_IN_ITERATIONS = 1000, MAX_ITERATIONS = 10000, DRAW_CYCLE = 20)
{
	# Author:Li Qinglong
	# Input:
	# y 	:an n*k matrix, observed data, each row is an individual's rank of items
	# p 	:number of parameters in MVNOS model
	# Z 	:a n*k*p array of covariates associated with all judges
	# beta0 :a 1*p matrix, prior normal distribution mean parameters
	# A0 	:a p*p matrix, prior normal distribution variance-covariance matrix
	# alpha :singular, prior Wishart distribution degree of freedom
	# P 	:a (k-1)*(k-1) matrix, prior Wishart distribution scale matrix

	# Output:
	# A list of Gibbs sampling trace

	# 	require(rjags)
	# Initialization
	# Prior distribution parameters
	# Default ones
	item_name = colnames(y)
	y = as.matrix(y)
	names(y) = NULL
	n = dim(y)[1] # number of individuls
	k = dim(y)[2] # number of items
	if (is.null(beta0)) beta0 = rep(0, p)
	
	if (is.null(A0)) A0 = diag(100, ncol = p, nrow = p)
	if (any(dim(A0) != c(p, p))) message("A0 shouble a p * p matrix.")
	
	if (is.null(alpha)) alpha = k + 1 
	
	if (is.null(P)) P = diag(k + 1, ncol = k - 1, nrow = k - 1)
	if (any(dim(P) != c(k - 1, k - 1))) message("P shouble a (k-1) * (k-1) matrix.")

	X = array(0, dim = c(n, k - 1, p))
	for (j in 1:n)
	{
	    Zk = Z[j, k,]
	    X[j, , ] = t(t(Z[j, 1:(k-1), ]) - Z[j, k, ])
	}
	
	#Starting value of w(Standardized rank score)
	yk = matrix(rep(y[, k], k), ncol = k)	
	w = (y - yk) / sqrt((k ^ 2 - 1) / 12)

	##################################################################
	# Start to using JAGS
	# Total iterates of Gibbs sampling
	MAX_ITERATIONS = MAX_ITERATIONS
	# Number of iterates to be reduced(burnt in)
	BURN_IN_ITERATIONS = BURN_IN_ITERATIONS
	data <- list(X = X, y = y, n = n, p = p, k = k, alpha = alpha, beta0 = beta0, A0 = A0, P = P)
	init <- list(w = w)
	init$w[, k] = NA

	# JAGS code
	# When p=1, beta ~ dnorm(beta0, A0)
	if (p == 1) 
	{
		strModelCode = "
		var X[n, k - 1, p], bounds[n, k - 1, 2], beta[p]

		data
		{
			for(i in 1:n)
			{
				for(j in 1:(k-1))
				{
					ones[i, j] <- 1
				}
			}
			lower <- -1e+5
			upper <-  1e+5
		}

		model
		{
			for (i in 1:n)
			{
				for (j in 1:(k-1))
				{
					bounds[i, j, 1] <- equals(y[i,j], 1) * lower + inprod(w[i, ], equals(y[i, ], y[i, j] - 1))
					bounds[i, j, 2] <- equals(y[i,j], k) * upper + inprod(w[i, ], equals(y[i, ], y[i, j] + 1))
		         	ones[i, j] ~ dinterval(w[i, j], bounds[i, j, ])
		        }
		        w[i, 1:(k-1)] ~ dmnorm(X[i, , ] * beta, G)
		      	w[i, k] <- 0
			}
			beta ~ dnorm(beta0, A0)
		   	G ~ dwish(P, alpha)
			Sigma <- inverse(G)
		}" 
	} else
	{
		strModelCode = "
		var X[n, k - 1, p], bounds[n, k - 1, 2], beta[p]

		data
		{
			for(i in 1:n)
			{
				for(j in 1:(k-1))
				{
					ones[i, j] <- 1
				}
			}
			lower <- -1e+5
			upper <-  1e+5
		}

		model
		{
			for (i in 1:n)
			{
				for (j in 1:(k-1))
				{
					bounds[i, j, 1] <- equals(y[i,j], 1) * lower + inprod(w[i, ], equals(y[i, ], y[i, j] - 1))
					bounds[i, j, 2] <- equals(y[i,j], k) * upper + inprod(w[i, ], equals(y[i, ], y[i, j] + 1))
		         	ones[i, j] ~ dinterval(w[i, j], bounds[i, j, ])
		        }
		        w[i, 1:(k-1)] ~ dmnorm(X[i, , ] %*% beta, G)
		      	w[i, k] <- 0
			}
			beta ~ dmnorm(beta0, A0)
		   	G ~ dwish(P, alpha)
			Sigma <- inverse(G)
		}"
	}
	jags <- rjags::jags.model(textConnection(strModelCode), data, init)
	update(jags, BURN_IN_ITERATIONS)
	samp <- rjags::coda.samples(jags, c("beta", "Sigma"), MAX_ITERATIONS)
	# End of running JAGS
	#################################################################

	posterior_data = as.matrix(samp)
	draw_index = seq(DRAW_CYCLE, MAX_ITERATIONS, DRAW_CYCLE)
	posterior_trace = posterior_data[draw_index, ]

	Sigma_trace = posterior_trace[, 1:(k-1)^2]
	beta_trace = posterior_trace[, ((k - 1)^2 + 1):((k - 1)^2 + p)]
	# Scaled by Sigma11
	Sigma11 = Sigma_trace[, 1]
	beta_trace = as.matrix(beta_trace / sqrt(Sigma11))
	Sigma_trace = Sigma_trace / Sigma11
	# Create an output list
	output_list = summary_trace(beta_trace, Sigma_trace,  item_name = item_name)
	return(output_list)
}