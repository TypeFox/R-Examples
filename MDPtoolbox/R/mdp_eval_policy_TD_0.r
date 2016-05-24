# function [V,mean_discrepancy] = mdp_eval_policy_TD_0(P,R,discount,policy,N)

mdp_eval_policy_TD_0 <- function(P,R,discount,policy,N) {

#  mdp_eval_policy_TD_0  Evaluation of the value function, using the TD(0) algorithm 
# 
#  Arguments
#  -------------------------------------------------------------------------
#  Let S = number of states, A = number of actions
#    P(SxSxA)  = transition matrix 
#               P could be an array with 3 dimensions or 
#               a cell array (1xA), each cell containing a matrix (SxS) possibly sparse
#    R(SxSxA) or (SxA) = reward matrix
#               R could be an array with 3 dimensions (SxSxA) or 
#               a cell array (1xA), each cell containing a sparse matrix (SxS) or
#               a 2D array(SxA) possibly sparse  
#    discount  = discount rate in ]0; 1[
#    policy(S) = optimal policy
#    N(optional) = number of iterations to execute, default value: 10000.
#               It is an integer greater than the default value. 
#  Evaluation --------------------------------------------------------------
#    V(S)   = optimal value function.
#    mean_discrepancy(N/100) = vector of V discrepancy mean over 100 iterations
#              Then the length of this vector for the default value of N is 100.

# check of arguments

if ( nargs() < 4 | nargs() > 5 ) {
	print('--------------------------------------------------------')
	print('The number of arguments must be 4 or 5')
	print('--------------------------------------------------------')
} else if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( nargs() == 5 & ifelse(!missing(N), N < 10000, F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: N must be upper than 10000')
	print('--------------------------------------------------------')
} else {
	# initialization of optional argument
	
	if (nargs() < 5) {
		N <- 10000
	}
	
	# Assigning the number of states and actions
	if (is.list(P)) {
		S <- dim(P[[1]])[1]
		A <- length(P)
	} else {
		S <- dim(P)[1]
		A <- dim(P)[3]
	}
	
	# Initializations
	s <- sample(1:S, 1, replace=T)		# Initial state choice
	V <- numeric(S)
	mean_discrepancy <- NULL		# vector of mean V variations
	discrepancy <- NULL			# vector of V variation
	
	for (n in 1:N) {
		# Reinitialisation of trajectories every 100 transitions
		if ( n %% 100 == 0 ) {
			s <- sample(1:S, 1, replace=T)
		}
		
		# Select an action: here action of the policy
		a <- policy[s]
		
		# Simulate next state s_new
		p_s_new <- runif(1)
		p <- 0
		s_new <- 0
		
		while ((p <= p_s_new) & (s_new < S)) {
			s_new <- s_new + 1
			if (is.list(P)) {
				p <- p + P[[a]][s,s_new]
			} else {
				p <- p + P[s,s_new,a]
			}
		}
		
		# Update V
		if (is.list(R)) {
			r <- R[[a]][s,s_new]
		} else {
			if(length(dim(R)) == 3) {
				r <- R[s,s_new,a]
			} else {
				r <- R[s,a]
			}
		}
		
		delta <- r + discount*V[s_new] - V[s]
		dV <- (1/sqrt(n+2))*delta
		V[s] <- V[s] + dV
		
		# Update current state
		s <- s_new
		
		# Memorize mean V variations on each trajectory (100 transitions)
		discrepancy[(n %% 100) + 1] = abs(dV)
		if(length(discrepancy) == 100) {
			mean_discrepancy <- c(mean_discrepancy, mean(discrepancy))
			discrepancy <- NULL
		}
	}
	return(list("V"=V,"mean_discrepancy"=mean_discrepancy))
}

}
