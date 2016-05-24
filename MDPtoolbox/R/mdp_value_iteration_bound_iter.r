mdp_value_iteration_bound_iter <- function(P, R, discount, epsilon, V0) {

# check of arguments
if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( nargs() > 3 & ifelse(!missing(epsilon), ifelse(epsilon < 0, T, F), F)  ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: epsilon must be upper than 0')
	print('--------------------------------------------------------')
} else if ( is.list(P) & nargs() > 4 & ifelse(!missing(V0), ifelse(length(V0) != dim(P[[1]])[1], T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: V0 must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( !is.list(P) & nargs() > 4 & ifelse(!missing(V0), ifelse(length(V0) != dim(P)[1], T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: V0 must have the same dimension as P')
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
	
	# set default values
	if (nargs() < 5) {
		V0 <- numeric(S)
	}
	if (nargs() < 4) {
		epsilon <- 0.01
	}

# See Markov Decision Processes, M. L. Puterman, 
# Wiley-Interscience Publication, 1994 
# p 202, Theorem 6.6.6
# k =    max     [1 - S min[ P(j|s,a), p(j|s',a')] ]
#     s,a,s',a'       j

	k <- 0
	h <- numeric(S)
	
	if (is.list(P)) {
		for (ss in 1:S) {
			PP <- NULL
			for (tt in 1:A) {
				PP <- rbind(PP,P[[tt]][,ss])
			}
			h[ss] <- min(min(PP))
		}
	} else {
		for (ss in 1:S) {
			h[ss] <- min(min(P[,ss,]))
		}
	}
	k <- 1 - sum(h)
	V1 <- mdp_bellman_operator(P,PR,discount,V0)[[1]]
	# p 201, Proposition 6.6.5
	max_iter <- log ( (epsilon*(1-discount)/discount) / mdp_span(V1-V0) ) / log(discount*k)

}

max_iter <- ceiling(max_iter)

return(max_iter)

}
