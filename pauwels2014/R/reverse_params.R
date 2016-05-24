if(getRversion() >= "2.15.1")  utils::globalVariables("knobjs")
reverse_params <-
function(true_params, params, transform_params){
	## Find theta such that transform_params(theta) = true_params
	## This is done by function maximization
	
	knobj <- c()

	temp <- rep(50, length(params))
	names(temp) <- names(params)
	knobj$global_parameters$tol <- 1e-20
	knobj$global_parameters$max_it = 200
	knobj$global_parameters$beta = 2
	knobj$global_parameters$c = 0.0001
	knobj$transform_params <- transform_params


	BFGS_special(
		temp, 
		fun_like = function(x, knobj, fail_incoming = F){
		print(x)
			temp <- -sum((knobj$transform_params(x)- true_params)^2); 
			if(fail_incoming){
				res <- c(); 
				res$res <- temp; 
				res$fail <- F
			}else{
				res <- temp
			}
			res
		}, 
		knobjs[[1]], verbose = T
	)$theta
}
