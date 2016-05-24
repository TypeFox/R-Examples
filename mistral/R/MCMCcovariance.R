## -----------------------------------------------------------------------------
## Fonction MCMCcovariance
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

MCMCcovariance = function(samples, n_seeds, chain_length, VA_function, VA_values=apply(samples,2,VA_function), N=length(VA_values), VA_esp, VA_var){

	if(missing(VA_esp)) VA_esp = mean(VA_values)

	#the i_th element of the j_th chain is located at (i-1)*n_seeds+j in G
	if(length(VA_values)<(chain_length*n_seeds)) {VA_values[(length(VA_values)+1):(chain_length*n_seeds)] = NA}
	dim(VA_values) = c(chain_length,n_seeds)
	VA_values = t(VA_values)

	R_tmp = NA*(1:chain_length)
	shift_max = chain_length-1
	for (shift in 0:shift_max) {
		k = 1
		s = 0
		n_sum = 0
		while((k+shift)<=chain_length) {
			s = sum(VA_values[,k]*VA_values[,k+shift],na.rm=TRUE) + s
			n_sum = n_sum + sum(!is.na(VA_values[,k]*VA_values[,k+shift]))
			k = k+1
		}
		R_tmp[shift+1] = s/n_sum - VA_esp^2
	}
	if(missing(VA_var)){VA_var = R_tmp[1]; cat(" VA_var =",VA_var,"\n")}# = R0
	R = R_tmp[-1]

	cat("#Calculate gamma \n")
	MC_gamma = 2*sum((seq(f=1,t=1,l=shift_max)-c(1:shift_max)*n_seeds/N)*R/VA_var)
	cat(" MC_gamma =",MC_gamma,"\n")

	cat("#Calculate Monte-Carlo variance\n")
	MC_var = 1/N*(VA_var+2*sum((seq(f=1,t=1,l=shift_max)-c(1:shift_max)*n_seeds/N)*R))
	cat(" MC_var =",MC_var,"\n")

	cat("#Calculate Monte-Carlo cov\n")
	MC_delta = sqrt(MC_var)/VA_esp
	cat(" MC_delta =",MC_delta,"\n")

	res = list(gamma=MC_gamma,var=MC_var,cov=MC_delta)
} 
