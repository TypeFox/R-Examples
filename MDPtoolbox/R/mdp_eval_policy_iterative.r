# function Vpolicy = mdp_eval_policy_iterative(P, R, discount, policy, V0, epsilon, max_iter)

mdp_eval_policy_iterative <- function(P, R, discount, policy, V0, epsilon, max_iter) {

#  mdp_eval_policy_iterative   Policy evaluation using iteration. 
#  Arguments -------------------------------------------------------------
#  Let S = number of states, A = number of actions
#    P(SxSxA)  = transition matrix 
#                P could be an array with 3 dimensions or 
#                a cell array (1xS), each cell containing a matrix possibly sparse
#    R(SxSxA) or (SxA) = reward matrix
#               R could be an array with 3 dimensions (SxSxA) or 
#               a cell array (1xA), each cell containing a sparse matrix (SxS) or
#               a 2D array(SxA) possibly sparse  
#    discount  = discount rate in ]0; 1[
#    policy(S) = a policy
#    V0(S)     = starting value function, optional (default : zeros(S,1))
#    epsilon   = epsilon-optimal policy search, upper than 0,
#                optional (default : 0.0001)
#    max_iter  = maximum number of iteration to be done, upper than 0, 
#                optional (default : 10000)
#  Evaluation -------------------------------------------------------------
#    Vpolicy(S) = value function, associated to a specific policy
# --------------------------------------------------------------------------
#  In verbose mode, at each iteration, displays the condition which stopped iterations:
#  epsilon-optimum value function found or maximum number of iterations reached.

# VERBOSE STUFF

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
} else if ( is.list(P) & nargs() > 4 & ifelse(!missing(V0), ifelse(length(V0) != dim(P[[1]])[1], T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: V0 must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( !is.list(P) & nargs() > 4 & ifelse(!missing(V0), ifelse(length(V0) != dim(P)[1], T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: V0 must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( nargs() > 5 & ifelse(!missing(epsilon), ifelse(epsilon < 0, T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: epsilon must be upper than 0')
	print('--------------------------------------------------------')
} else if ( nargs() > 6 & ifelse(!missing(max_iter), ifelse(max_iter <= 0, T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: The maximum number of iteration must be upper than 0')
	print('--------------------------------------------------------')
} else {
	# initialization of optional arguments
	if (is.list(P)) {
		S <- dim(P[[1]])[1]
	} else {
		S <- dim(P)[1]
	}
	if (nargs() < 5) {
		V0 <- numeric(S)
	}
	if (nargs() < 6) {
		epsilon <- 0.0001
	}
	if (nargs() < 7) {
		max_iter <- 10000
	}
	
	compute <- mdp_computePpolicyPRpolicy(P, R, policy)
	Ppolicy <- compute[[1]]
	PRpolicy <- compute[[2]]
	
	iter <- 0
	Vpolicy <- V0
	is_done <- F
	
	while (!is_done) {
		iter <- iter + 1
		Vprev <- Vpolicy
		Vpolicy <- PRpolicy + discount * Ppolicy %*% Vprev
		variation <- max(abs(Vpolicy - Vprev))
		if (variation < epsilon) {
			is_done <- T
			# VERBOSE STUFF
			print('MDP Toolbox: iterations stopped, epsilon-optimal value function')
		} else if (iter == max_iter) {
			is_done <- T
			print(paste('MDP Toolbox: iterations stopped by maximum number of iteration condition', max_iter))
		}
	}
}

return(Vpolicy)

}
