mdp_relative_value_iteration <- function(P, R, epsilon, max_iter) {

start<-as.POSIXlt(Sys.time())

# VERBOSE STUFF

# check of arguments

if ( nargs() > 3 & ifelse(!missing(epsilon), epsilon < 0, F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: epsilon must be upper than 0')
	print('--------------------------------------------------------')
} else if ( nargs() > 4 & ifelse(!missing(max_iter), max_iter <= 0, F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: The maximum number of iteration must be upper than 0')
	print('--------------------------------------------------------')
} else {

	if (is.list(P)) {
		S <- dim(P[[1]])[1]
		A <- length(P)
	} else {
		S <- dim(P)[1]
		A <- dim(P)[3]
	}
	
	PR <- mdp_computePR(P,R)
	
	# set default values
	if (nargs() < 4) {
		max_iter <- 1000
	}
	if (nargs() < 3) {
		epsilon <- 0.01
	}
	
	U <- numeric(S)
	
	# VERBOSE STUFF
	
	iter <- 0
	is_done <- F
	Q1 <- numeric(A)
	Q2 <- matrix(0, S, A)
	
	while (!is_done) {
		iter <- iter + 1
		
		if (is.list(P)) {
			for ( a in 1:A) {
				Q1[a] <- PR[S,a] + P[[a]][S,]%*%U
			}
		} else {
			for ( a in 1:A) {
				Q1[a] <- PR[S,a] + P[S,,a]%*%U
			}
		}
		
		g <- max(Q1)
		
		if (is.list(P)) {
			for ( a in 1:A) {
				Q2[,a] <- as.matrix(PR[,a] + P[[a]]%*%U)
			}
		} else {
			for ( a in 1:A) {
				Q2[,a] <- PR[,a] + P[,,a]%*%U
			}
		}
		
		Unext <- apply(Q2, 1, max)
		policy <- apply(Q2, 1, which.max)
		
		Unext <- Unext - g
		variation <- mdp_span(Unext-U)
		
		# VERBOSE STUFF
		
		if (variation < epsilon) {
			is_done <- T
			
			# VERBOSE STUFF
			print('MDP Toolbox: iterations stopped, epsilon-optimal policy found')
		} else if (iter == max_iter) {
			is_done <- T
			print('MDP Toolbox: iterations stopped by maximum number of iteration condition')
		} else {
			U <- Unext
		}
	}
}

end <-as.POSIXlt(Sys.time())

return(list(U, policy, g, end-start))

}

