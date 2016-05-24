WeibullIntegrate <- function(x, x_i_noncens=NULL, density, param_y_i, param_delta_i, param_lambda, param_gamma, param_beta, intlimit = 10^-10, ForIntegrate=TRUE){
	if(is.null(x_i_noncens)==FALSE){
		x_i <- matrix(nrow=1, ncol=(length(x_i_noncens)+1))
		for(ii in 1:length(x_i_noncens)){
			x_i[ii] <- x_i_noncens[ii]
		}
		x_i[length(x_i)] <- x
	}
	if(is.null(x_i_noncens)==TRUE){
		x_i <- x
	}
    
	#result <- (f_conditional(t=param_y_i,z=x_i,lambda=param_lambda,a=param_gamma,beta=param_beta)^param_delta_i) *
	#		   (S_conditional(t=param_y_i,z=x_i,lambda=param_lambda,a=param_gamma,beta=param_beta)^(1-param_delta_i)) *
	#		   density(x)
	
	h_i_cond <- h_conditional(t=param_y_i,z=x_i,lambda=param_lambda,a=param_gamma,beta=param_beta)^param_delta_i
	S_i_cond <- S_conditional(t=param_y_i,z=x_i,lambda=param_lambda,a=param_gamma,beta=param_beta)
	if(h_i_cond==Inf & S_i_cond==0){
		result <- 0
	}
	if(h_i_cond==-Inf & S_i_cond==0){
		result <- 0
	}
	if(h_i_cond!=Inf | S_i_cond!=0){
    if(h_i_cond!=-Inf | S_i_cond!=0){
        result <- (h_conditional(t=param_y_i,z=x_i,lambda=param_lambda,a=param_gamma,beta=param_beta)^param_delta_i) *
                    S_conditional(t=param_y_i,z=x_i,lambda=param_lambda,a=param_gamma,beta=param_beta) *
                    density(x)
        
		if(ForIntegrate==TRUE){
			if(result <= intlimit){result <- 0}
		}
	}}
    return(result)
}















#