survtime <- function(param_ts, covariate, lambda, gamma, beta){
     
	result_survtime <- exp((-exp(covariate%*%beta)) * lambda * (param_ts^gamma))
     
     return(result_survtime)
}











#