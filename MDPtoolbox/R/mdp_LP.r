mdp_LP <- function(P, R, discount) {

start<-as.POSIXlt(Sys.time())

if ( discount <= 0 | discount > 1 ) {
	print('--------------------------------------------------------')
	print('MDP Toolbox ERROR: Discount rate must be in ]0; 1]')
	print('--------------------------------------------------------')
} else if (is.list(P)) {
	print('--------------------------------------------------------')
	print('MDPR ERROR: mdp_LP cannot be used with sparse matrices')
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

# The objective is to resolve : min V such as V >= PR + discount*P*V
# The function linprog of the optimisation Toolbox of Mathworks resolves :
# min f'*x such as M * x <= b
# So the objective could be expressed as : min V / (discount*P-I) * V <= - PR
# To avoid loop on states, the matrix M is structured following actions M(A*S,S)

	f <- rep(1,S)
	M <- NULL
	
	if (is.list(P)) {
		for (a in 1:A) {
			M <- rbind2(M, discount*P[[a]] - Diagonal(S))
		}
	} else {
		for (a in 1:A) {
			M <- rbind(M, as.matrix(discount*P[,,a] - Diagonal(S)))
		}
	}
	
	# MATLAB: linprog(f,A,b) solves min f'*x such that A*x ≤ b.
	# R: solveLP(f,b,A) solves min f'*x, subject to A*x <= b and x >= 0.
  # Tricky with the number of rows of b and the size of A
	V <- as.numeric(solveLP(f,as.vector(-PR),M)$solution)
	
	bellman <- mdp_bellman_operator(P,PR,discount,V)
	V <- bellman[[1]]
	policy <- bellman[[2]]
	
	end <-as.POSIXlt(Sys.time())
	return(list("V"=V, "policy"=policy, "time"=end-start))

}

}
