armijo <-
function(theta, fun, p, gradf, valf, l = 1, beta = 2, c = 0.0001, ...){
	## Implements armijo condition
	
	temp <- sum( gradf * p )
	
	theta1 <- theta + l * p
	names(theta1) <- names(theta)
	valf1 <- fun(theta1, ...)
	temp1 <- (valf1 >= valf + c * l * temp) 
	if (is.na(temp1)){
		temp1 <- FALSE
	}
	
	if(temp1){
		next_it <- T
		valf11 <- valf
		while(next_it){
			l_prev <- l
			l <- l*beta
			theta1 <- theta + l * p
			names(theta1) <- names(theta)
			valf1 <- valf11
			valf11 <- fun(theta1, ...)
			next_it <- ((valf11 >= valf + c * l * temp) & (l>1e-20))
			if(is.na(next_it)){
				next_it <- FALSE
			}
		}
		l <- l_prev
	}else{
		next_it <- T
		while(next_it){
			l <- l / beta
			theta1 <- theta + l * p
			names(theta1) <- names(theta)
			valf1 <- fun(theta1, ...)
			next_it <- ((valf1 < valf + c * l * temp) & (l>1e-20))
			if(is.na(next_it)){
				next_it <- FALSE
			}
		}
	}
	l
}
