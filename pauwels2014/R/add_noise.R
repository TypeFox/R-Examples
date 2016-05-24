add_noise <-
function(data_theta_Ts){
	## Add heteroscedastic noise to simulations
	temp <- cbind(data_theta_Ts[,1], matrix(rnorm(n = (ncol(data_theta_Ts) - 1) * nrow(data_theta_Ts), mean = as.numeric(data_theta_Ts[,-1]), sd = sqrt(0.01 + 0.04 * as.numeric(data_theta_Ts[,-1])^2 ) ), nrow(data_theta_Ts), ncol(data_theta_Ts) - 1 ) )
	dimnames(temp) <- dimnames(data_theta_Ts)
	temp
}
