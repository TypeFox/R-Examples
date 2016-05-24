mdp_finite_horizon <- function(P, R, discount, N, h) {

start<-as.POSIXlt(Sys.time())

# VERBOSE STUFF

# check of arguments
if (N < 1) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: N must be upper than 0')
	print('--------------------------------------------------------')
} else if (discount <= 0 | discount > 1) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
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
	V <- matrix(0,S,N+1)
	policy <- matrix(0,S,N)
	
	if (nargs() == 5) {
	#	V[,N+1] <- h
	}
	
	PR <- mdp_computePR(P,R)
	
	for (n in 0:(N-1)) {
		bellman <- mdp_bellman_operator(P,PR,discount,V[,N-n+1])
		W <- bellman[[1]]
		X <- bellman[[2]]
		V[,N-n] <- W
		policy[,N-n] <- X
	}
}

end <-as.POSIXlt(Sys.time())

return(list("V"=V, "policy"=policy, "time" = end-start))

}
