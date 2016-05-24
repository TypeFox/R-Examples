# function Vpolicy = mdp_eval_policy_matrix(P, R, discount, policy)

mdp_eval_policy_matrix <- function(P, R, discount, policy) {

# check of arguments

if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( ifelse(is.list(P), length(policy) != dim(P[[1]])[1], F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: policy must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( ifelse(!is.list(P), length(policy) != dim(P)[1], F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: policy must have the same dimension as P')
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
	
	compute <- mdp_computePpolicyPRpolicy(P, R, policy)
	Ppolicy <- compute[[1]]
	PRpolicy <- compute[[2]]
	
	# A comment:
	# V = PR + gPV  => (I-gP)V = PR  => V = inv(I-gP)*PR
	
	Vpolicy <- solve((Diagonal(S) - discount*Ppolicy), PRpolicy) # Check the left division
	
	return(as.numeric(Vpolicy))

}

}
