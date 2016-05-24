# Created by Matthew A. Birk
# Calculates pooled summary statistics
# Last updated: Jun 2015

#' Pooled Summary Descriptive Statistics
#'
#' Pools summary statistics when given mean and (optionally) a measurement of variability (choose one among \code{var}, \code{sd}, and \code{se}).
#'
#' @param mean numeric. A vector of mean values to be pooled.
#' @param n numeric. A vector of n values to be pooled.
#' @param var numeric. A vector of variance values to be pooled.
#' @param sd numeric. A vector of standard deviation values to be pooled.
#' @param se numeric. A vector of standard error of the mean vlaues to be pooled.
#'
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#' @seealso \code{\link{weighted.mean}}, \code{\link{se}}
#'
#' @examples
#' summ_stat(mean = c(0.68, 0.67), n = c(4, 5), sd = c(0.11, 0.15))
#' summ_stat(mean = 0.68, n = 3, se = 5)
#' summ_stat(mean = rnorm(1e4), n = rep(1, 1e4)) # Find pooled mean when variability is unknown.
#' 
#' @encoding UTF-8
#' @export
#' @import stats

summ_stat = function(mean, n, var, sd, se){
	pooled_mean = stats::weighted.mean(mean, n)
	if(missing(var)){
		if(missing(sd)){
			if(missing(se)) return(list(pooled_mean = pooled_mean))
			sd = se * sqrt(n)
		}
		var = sd ^ 2
	}
	N = sum(n)
	pooled_var = 1 / (N - 1) * (sum((n - 1) * var) + sum(n * (mean - weighted.mean(mean, n)) ^ 2))
	pooled_sd = sqrt(pooled_var)
	pooled_se = pooled_sd / sqrt(N)
	return(list(pooled_mean = pooled_mean, pooled_var = pooled_var, pooled_sd = pooled_sd, pooled_se = pooled_se))
}