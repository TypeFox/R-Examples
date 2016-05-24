LoglikWeibullSurvRegCens <- function(x, data_y, data_delta_loglik, data_cov_noncens=NULL, data_cov_cens, density, data_r_loglik, data_lowerbound, intlimit = 10^-10){
    
    n <- length(data_y)
    
    beta <- matrix(nrow = (length(x) - 2), ncol = 1)
    for(ii in 1:(length(x)-2)){
        beta[ii] <- x[ii+2]
    }
    
    lambda <- x[1]
    gamma <- x[2]
    
    l_beta <- 0
    
    if(lambda <= 0 | gamma <= 0){ # Because the parameters lambda and gamma cannot be negative!
        l_beta <- -Inf
    }
    
    if(lambda > 0 & gamma >0){
        for(ll in 1:n){
            r_i <- data_r_loglik[ll]
            delta_i <- data_delta_loglik[ll]
            if(is.null(data_cov_noncens)==TRUE){
                x_i_noncens <- NULL
            }
            if(is.null(data_cov_noncens)==FALSE){
                if(ncol(as.matrix(data_cov_noncens))==1){
                    x_i_noncens <- data_cov_noncens[ll]
                }
                if(ncol(as.matrix(data_cov_noncens))>1){
                    x_i_noncens <- data_cov_noncens[ll,]
                }
            }
            x_i_cens <- data_cov_cens[ll]
            y_i <- data_y[ll]
            if(r_i == 1){
                # The first two "if" are just for the case where there are no censored observation, because then "data_cov_cens[data_r_loglikelihood==0 & data_cov_cens<=x_i_cens]" gives you an empty answer and you get "warnings" !
                if(length(data_cov_cens[data_r_loglik==0 & data_cov_cens<=x_i_cens])==0){
                    c <- -Inf
                }
                if(length(data_cov_cens[data_r_loglik==0 & data_cov_cens<=x_i_cens])!=0){
                    c <- max(data_cov_cens[data_r_loglik==0 & data_cov_cens<=x_i_cens])
                }
                pi_i <- CDF(c, density=density)
                if(pi_i >= 1 & pi_i <= 1.00001){pi_i <- 1} # Because there are some approximations errors when computing the integration of the density function
                f_R_i_1 <- dbinom(1, size=1, prob=pi_i)
                
                l_beta <- l_beta + log(f_R_i_1 * WeibullIntegrate(x=x_i_cens, x_i_noncens=x_i_noncens, density=density, param_y_i=y_i, param_delta_i=delta_i,
                param_lambda=lambda, param_gamma=gamma, param_beta=beta, intlimit=intlimit, ForIntegrate=FALSE))
            }
            if(r_i == 0){
                pi_i <- CDF(x_i_cens, density=density)
                if(pi_i >= 1 & pi_i <= 1.00001){pi_i <- 1} # Because there are some approximations errors when computing the integration of the density function
                f_R_i_0 <- dbinom(0, size=1, prob=pi_i)
                if(is.na(data_lowerbound[ll])==FALSE){
                    LowerBound <- data_lowerbound[ll]
                }
                if(is.na(data_lowerbound[ll])==TRUE){
                    LowerBound <- -Inf
                }
                beta_noncens <- beta[1:(length(beta)-1)]
                l_beta <- l_beta + log(integrate(f=Vectorize(WeibullIntegrate, vectorize.args=c("x")), lower=LowerBound, upper=x_i_cens, x_i_noncens=x_i_noncens, density=density,
                param_delta_i=delta_i, param_beta=beta, param_lambda=lambda, param_gamma=gamma, param_y_i=y_i, intlimit=intlimit)$value *
                f_R_i_0)
            }
        }
    }
    
    if(l_beta == -Inf){l_beta <- -100000}
    
    return(l_beta)
}













#