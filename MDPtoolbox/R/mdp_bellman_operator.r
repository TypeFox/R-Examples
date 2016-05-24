mdp_bellman_operator <- function(P, PR, discount, Vprev) {

if (discount <= 0 | discount > 1) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if ( is.list(P) & ifelse(!missing(Vprev), ifelse(length(Vprev) != dim(P[[1]])[1], T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Vprev must have the same dimension as P')
	print('--------------------------------------------------------')
} else if ( !is.list(P) & ifelse(!missing(Vprev), ifelse(length(Vprev) != dim(P)[1], T, F), F) ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Vprev must have the same dimension as P')
	print('--------------------------------------------------------') 
} else {
	Q <- matrix(0, dim(PR)[1], dim(PR)[2])
	if (is.list(P)) {
		A <- length(P)
		for (a in 1:A) {
    # Here loose sparse property - but I guess it must be done anyhow, note that as.matrix is called on the results, not on P[a]
			Q[,a] <- as.matrix(PR[,a] + discount*P[a][[1]] %*% Vprev)
		}
	} else {
		A <- dim(P)[3]
		for (a in 1:A) {
			Q[,a] <- PR[,a] + discount*P[,,a] %*% Vprev
		}
	}
	return(list( "V"=apply(Q, 1, max), "policy"=apply(Q, 1, which.max) ))
}

}
