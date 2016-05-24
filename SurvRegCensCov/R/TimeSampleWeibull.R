TimeSampleWeibull <- function(covariate_noncens = NULL, covariate_cens, lambda, gamma, beta){
     
  covariate_noncens <- as.matrix(covariate_noncens)
  
    if(is.null(covariate_noncens) == FALSE){covariate <- cbind(covariate_noncens, covariate_cens)}
    if(is.null(covariate_noncens) == TRUE){covariate <- covariate_cens}
    
	U_simulate <- runif(length(covariate_cens))
	to_solve <- function(x, covariate=covariate, lambda, gamma, beta, u_unif){
		result_to_solve <- F_time(param_ts=x, covariate=covariate, lambda=lambda, gamma=gamma, beta=beta) - u_unif
        return(result_to_solve)
	}
	result_generate_time_sample_weibull <- matrix(ncol=1, nrow=length(covariate_cens))
	for(ii in 1:length(covariate_cens)){
        if(is.null(covariate_noncens) == FALSE){
            result_generate_time_sample_weibull[ii] <- uniroot(f=to_solve, interval=c(0, 100000), covariate=covariate[ii,], lambda=lambda, gamma=gamma, beta=beta, u_unif=U_simulate[ii])$root
        }
        if(is.null(covariate_noncens) == TRUE){
            result_generate_time_sample_weibull[ii] <- uniroot(f=to_solve, interval=c(0, 100000), covariate=covariate[ii], lambda=lambda, gamma=gamma, beta=beta, u_unif=U_simulate[ii])$root
        }
	}
    return(result_generate_time_sample_weibull)
}











#