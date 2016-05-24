F_time <- function(param_ts, covariate, lambda, beta, gamma){
	result_F_time <- 1 - survtime(param_ts = param_ts, covariate = covariate, lambda = lambda, gamma = gamma, beta = beta)
    return(result_F_time)
}











#