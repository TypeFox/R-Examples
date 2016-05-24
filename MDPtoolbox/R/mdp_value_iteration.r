mdp_value_iteration <- function(P, R, discount, epsilon, max_iter, V0) {

start<-as.POSIXlt(Sys.time())

# VERBOSE STUFF

# check of arguments

if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( nargs() > 3 & ifelse(!missing(epsilon), ifelse(epsilon < 0, T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: epsilon must be upper than 0')
	print('--------------------------------------------------------')
} else if ( nargs() > 4 & ifelse(!missing(max_iter), ifelse(max_iter <= 0, T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: The maximum number of iteration must be upper than 0')
	print('--------------------------------------------------------')
} else if ( is.list(P) & nargs() > 5 & ifelse(!missing(V0), ifelse(length(V0) != dim(P[[1]])[1], T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: V0 must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( !is.list(P) & nargs() > 5 & ifelse(!missing(V0), ifelse(length(V0) != dim(P)[1], T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: V0 must have the same dimension as P')
	print('--------------------------------------------------------')
} else {
	if (discount == 1) {
		print('--------------------------------------------------------')
		print('MDP Toolbox WARNING: check conditions of convergence.')
		print('With no discount, convergence is not always assumed.')
		print('--------------------------------------------------------')
	}
	
	# Beware global scope!
	if (is.list(P)) {
		S <- dim(P[[1]])[1]
		A <- length(P)
	} else {
		S <- dim(P)[1]
		A <- dim(P)[3]
	}
	
	PR <- mdp_computePR(P,R)
	
	# initialization of optional arguments
	if (nargs() < 6) {
		V0 <- numeric(S)
	}
	if (nargs() < 4) {
		epsilon <- 0.01
	}

	# compute a bound for the number of iterations
	if (discount != 1) computed_max_iter <- 5000 # mdp_value_iteration_bound_iter(P, R, discount, epsilon, V0)

	if (nargs() < 5) {
		if (discount != 1) {
			max_iter <- computed_max_iter
		} else {
			max_iter <- 5000
		}
	} else {
		if (discount != 1 & max_iter > computed_max_iter) {
			print(paste('MDP Toolbox WARNING: max_iter is bounded by ', computed_max_iter))
			max_iter <- computed_max_iter
		}
	}
	
	# computation of threshold of variation for V for an epsilon-optimal policy
	if (discount != 1) {
		thresh <- epsilon * (1-discount)/discount
	} else {
		thresh <- epsilon
	}
	
	# VERBOSE STUFF
	
	iter <- 0
	V <- V0
	is_done <- F
	
	while(!is_done) {
		iter <- iter + 1
		Vprev <- V
		
		bellman <- mdp_bellman_operator(P,PR,discount,V)
		V <- bellman[[1]]
		policy <- bellman[[2]]
		
		variation <- mdp_span(V - Vprev)
		
		# VERBOSE STUFF
		
		if (variation < thresh) {
			is_done <- T
			
			# VERBOSE STUFF
			print('MDP Toolbox: iterations stopped, epsilon-optimal policy found')
		} else if (iter == max_iter) {
			is_done <- T
			print('MDP Toolbox: iterations stopped by maximum number of iteration condition')
		}
	}
}

end <-as.POSIXlt(Sys.time())

return(list("V"=V, "policy"=policy, "iter"=iter, "time"=end-start, "epsilon" = epsilon, "discount" = discount))

}
