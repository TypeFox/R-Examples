mdp_policy_iteration_modified <- function(P, R, discount, epsilon, max_iter) {

start<-as.POSIXlt(Sys.time())

# VERBOSE STUFF

# check of arguments

if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( nargs() > 4 & ifelse(!missing(epsilon), ifelse(epsilon < 0, T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: epsilon must be upper than 0')
	print('--------------------------------------------------------')
} else if ( nargs() > 5 & ifelse(!missing(max_iter), ifelse(max_iter <= 0, T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: The maximum number of iteration must be upper than 0')
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
	if (nargs() < 5) {
		max_iter <- 10
	}
	if (nargs() < 4) {
		epsilon <- 0.01
	}
	
	# computation of threshold of variation for V for an epsilon-optimal policy
	if (discount != 1) {
		thresh = epsilon * (1-discount)/discount
	} else {
		thresh = epsilon
	}
	
	if(discount == 1) {
		V <- numeric(S)
	} else {
		V <- 1/(1-discount)*min(min(PR))*rep(1,S)
	}
	
	# VERBOSE STUFF
	
	iter <- 0
	is_done <- F
	
	while(!is_done) {
		iter <- iter + 1
		bellman <- mdp_bellman_operator(P,PR,discount,V)
		Vnext <- bellman[[1]]
		policy <- bellman[[2]]
		
		compute <- mdp_computePpolicyPRpolicy(P, PR, policy)
		Ppolicy <- compute[[1]]
		PRpolicy <- compute[[2]]
		
		variation <- mdp_span(Vnext - V)
		V <- Vnext
		
		if (variation < thresh) {
			is_done <- T
		} else {
			# STRANGE VERBOSE STUFF
      V <- mdp_eval_policy_iterative(P, PR, discount, policy, V, epsilon, max_iter)
		}
	}
	
}

end <-as.POSIXlt(Sys.time())

return(list("V"=V, "policy"=policy, "iter"=iter, "time"=end-start))

}
