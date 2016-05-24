mdp_Q_learning <- function(P, R, discount, N) {

# check of arguments

if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( nargs() >= 4 & ifelse(!missing(N), N < 10000, F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: N must be upper than 10000')
	print('--------------------------------------------------------')
} else {
	# initialization of optional argument
	
	if (nargs() < 4) {
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
	Q <- matrix(0,S,A)
	dQ <- matrix(0,S,A)
	mean_discrepancy <- NULL
	discrepancy <- NULL

	# Initial state choice
	state <- sample(1:S, 1, replace=T)
	
	for (n in 1:N) {
	
		# Reinitialisation of trajectories every 100 transitions
		if ( n %% 100 == 0 ) {
			state <- sample(1:S, 1, replace=T)
		}
		
		# Action choice : greedy with increasing probability
		# probability 1-(1/log(n+2)) can be changed
		
		pn <- runif(1)
		if ( pn < (1-(1/log(n+2))) ) {
			optimal_action <- max(Q[state,])
			a <- which.max(Q[state,])
		} else {
			a <- sample(1:A, 1, replace=T)
		}
		
		# Simulating next state s_new and reward associated to <s,s_new,a>
		p_s_new <- runif(1)
		p <- 0
		s_new <- 0
		
		while ((p < p_s_new) & (s_new < S)) {
			s_new <- s_new + 1
			if (is.list(P)) {
				p <- p + P[[a]][state,s_new]
			} else {
				p <- p + P[state,s_new,a]
			}
		}
		
		if (is.list(R)) {
			r <- R[[a]][state,s_new]
		} else {
			if(length(dim(R)) == 3) {
				r <- R[state,s_new,a]
			} else {
				r <- R[state,a]
			}
		}
		
		# Updating the value of Q   
		# Decaying update coefficient (1/sqrt(n+2)) can be changed
		delta <- r + discount*max(Q[s_new,]) - Q[state,a]
		dQ <- (1/sqrt(n+2))*delta
		Q[state,a] <- Q[state,a] + dQ
		
		# Current state is updated
		state <- s_new
		
		# Computing and saving maximal values of the Q variation
		discrepancy[(n %% 100) + 1] = abs(dQ)
		if(length(discrepancy) == 100) {
			mean_discrepancy <- c(mean_discrepancy, mean(discrepancy))
			discrepancy <- NULL
		}
		
	}
	
	# compute the value function and the policy
	V <- apply(Q, 1, max)
	policy <- apply(Q, 1, which.max)

}

return(list("Q"=Q, "V"=V, "policy"=policy, "mean_discrepancy"=mean_discrepancy))

}

