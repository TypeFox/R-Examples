mdp_example_rand <- function(S, A, is_sparse, mask) {

# arguments checking
if (S < 1 | A < 1) {
	print('--------------------------------------------------------')
	print('MDPR ERROR: Both the number of states S ')
	print('and the number of actions A must be upper than 1.')
	print('--------------------------------------------------------')
} else if (nargs() == 4 & ifelse(!missing(mask), nrow(mask) != S | ncol(mask) != S, F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: mask must be a SxS matrix') 
	print('--------------------------------------------------------')
} else {
	# initialization of optional arguments
	
	if (nargs() < 3) {
		is_sparse <- F
	}
	
	if (nargs() < 4) {
		mask <- matrix(1,S,S)
	}
	
	if (is_sparse) {
		# definition of transition matrix : square stochastic matrix
		P <- list()
		for (a in 1:A) {
			PP <- Matrix(mask * runif(S*S), sparse = T)
			for (s in 1:S) {
				PP[s,] <- PP[s,] / sum( PP[s,] )
			}
		P[a] <- PP
		}
		# definition of reward matrix (values between -1 and +1)
		R <- list()
		for (a in 1:A) {
			R[a] <- Matrix(mask * (2*runif(S*S) - matrix(1,S,S)))
		}
	} else {
		# definition of transition matrix : square stochastic matrix
		P <- array(0, c(S,S,A))
		for (a in 1:A) {
			P[,,a] <- mask * runif(S*S)
			for (s in 1:S) {
				P[s,,a] <- P[s,,a] / sum( P[s,,a] )
			}
		}
		# definition of reward matrix (values between -1 and +1)
		R <- array(NA, c(S,S,A))
		for (a in 1:A) {
			R[,,a] <- mask * (2*runif(S*S) - matrix(1,S,S))
		}
	}
	
	return(list("P"=P,"R"=R))
	
}

}

