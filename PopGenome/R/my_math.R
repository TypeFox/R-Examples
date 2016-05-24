#berend krebs master thesis
#math background module



fstcal <- function(i,j){
	return (exp(BBB$d_alpha[i]+BBB$d_beta[j])/(1+exp(BBB$d_alpha[i]+BBB$d_beta[j])))
}


## log_prior_alpha
log_prior_alpha <- function(alpha_){
	return(log(0.5*(1/(BBB$sd_prior_alpha*sqrt(2*pi)))*exp(-(alpha_-BBB$m1_prior_alpha)*(alpha_-BBB$m1_prior_alpha)/(2*BBB$sd_prior_alpha^2)) + 0.5*(1/   (BBB$sd_prior_alpha*sqrt(2*pi)))*exp(-(alpha_-BBB$m2_prior_alpha)*(alpha_-BBB$m2_prior_alpha)/(2*BBB$sd_prior_alpha^2))))
}


