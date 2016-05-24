risk_theta_vect <-
function(thetas_trans, n_params){
	## Expected risk based on a sample in matrix form
	
	risk_current <- 0
	for(i in 1:nrow(thetas_trans)){
		risk_current <- risk_current + sum( log(t(t(thetas_trans) / thetas_trans[i,]) )^2 ) / n_params / nrow(thetas_trans)^2
	}
	risk_current
}
