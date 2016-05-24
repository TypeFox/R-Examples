compute_gradient <-
function(theta, fun, step = 1e-6, ...){
	## Finite difference method to compute the gradient of fun
	## at theta
	#ll <- fun(theta)
	temp <- sapply(1:length(theta), FUN = compute_gradient_coordinate, theta = theta, fun = function(x){fun(x,...)}, step = step)
	names(temp) <- names(theta)
	temp
}
