mdp_eval_policy_optimality <- function(P, R, discount, Vpolicy) {

# check of arguments

if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( ifelse(is.list(P), length(Vpolicy) != dim(P[[1]])[1], FALSE) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: policy must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( ifelse(!is.list(P), length(Vpolicy) != dim(P)[1], FALSE) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: policy must have the same dimension as P')
	print('--------------------------------------------------------')
} else {
	# compute Q(SxA)
	PR <- mdp_computePR(P,R)
	
	if (is.list(P)) {
		S <- dim(P[[1]])[1]
		A <- length(P)
		Q <- matrix(0, S, A)
		for ( a in 1:A) {
			Q[,a] <- as.matrix(PR[,a] + discount*P[[a]]%*%Vpolicy)
		}
	} else {
		S <- dim(P)[1]
		A <- dim(P)[3]
		Q <- matrix(0, S, A)
		for ( a in 1:A) {
			Q[,a] <- PR[,a] + discount*P[,,a]%*%Vpolicy
		}
	}
	
	# search near optimal actions a for each state s, satisfaying
	# Q(s,a) - Q(s, a*) < epsilon
	# where a* is the optimal action for state s
	
	epsilon <- 0.01
	optimal_actions <- abs(Q - matrix(rep(Vpolicy,A), length(Vpolicy), A)) < epsilon
	
	if ( max(colSums(optimal_actions)) == 1 ) {
		is_multiple <- F
	} else {
		is_multiple <- T
	}
	return(list("multiple"=is_multiple, "optimal_actions"=optimal_actions))
}

}
