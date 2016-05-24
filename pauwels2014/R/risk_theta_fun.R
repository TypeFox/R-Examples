risk_theta_fun <-
function(theta1, theta2, n_params){
	## Risk function
	1 / n_params * sum(log(theta1/theta2)^2)
}
