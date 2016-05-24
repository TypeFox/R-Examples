mdp_policy_iteration <- function(P, R, discount, policy0, max_iter, eval_type) {

start<-as.POSIXlt(Sys.time())

# VERBOSE STUFF

# check of arguments

if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( nargs() > 3 & is.list(P) & ifelse(!missing(policy0), length(policy0) != dim(P[[1]])[1], F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: policy must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( nargs() > 3 & !is.list(P) & ifelse(!missing(policy0), length(policy0) != dim(P)[1], F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: policy must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( nargs() > 4 & ifelse(!missing(max_iter), max_iter <= 0, F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: The maximum number of iteration must be upper than 0')
	print('--------------------------------------------------------')
} else {

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
		eval_type <- 0
	}
	if (nargs() < 5) {
		max_iter <- 1000
	}
	if (nargs() < 4) {
		# initialization of policy: 
		# the one wich maximizes the expected immediate reward
		bellman <- mdp_bellman_operator(P,PR,discount,numeric(S))
		Vunused <- bellman[[1]]
		policy0 <- bellman[[2]]
	}
	
	iter <- 0
	policy <- policy0
	is_done <- F
	
	while (!is_done) {
		iter <- iter + 1
		
		if (eval_type == 0) {
			V <- mdp_eval_policy_matrix(P,PR,discount,policy)
		} else {
			V <- mdp_eval_policy_iterative(P,PR,discount,policy)
		}
		
		bellman <- mdp_bellman_operator(P,PR,discount,V)
		Vnext <- bellman[[1]]
		policy_next <- bellman[[2]]
		
		n_different <- sum(policy_next != policy)
		
		# VERBOSE STUFF
		
		if ( setequal(policy_next,policy) | iter == max_iter ) {
			is_done <- T
		} else {
			policy <- policy_next
		}
	}
	
	end <-as.POSIXlt(Sys.time())
	return(list("V"=V, "policy"=policy, "iter"=iter, "time"=end-start))
}

}
