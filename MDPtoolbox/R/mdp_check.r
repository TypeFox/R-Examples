mdp_check <- function(P, R) {

is_error_detected <- FALSE
error_msg <- ''

# Check of P
if (is.list(P)) {
	s1 <- dim(P[[1]])[1]
	s2 <- dim(P[[1]])[2]
	a1 <- length(P)
} else {
	s1 <- dim(P)[1]
	s2 <- dim(P)[2]
	a1 <- dim(P)[3]
}

if (s1 < 1 | a1 < 1 | s1 != s2) {
	error_msg <- 'MDP Toolbox ERROR: The transition matrix must be on the form P(S,S,A) with S : number of states greater than 0 and A : number of action greater than 0'
	is_error_detected <- T
}

if (!is_error_detected) {
	a <- 1
	while (a <= a1) {
		if (is.list(P)) {
			 error_msg <- mdp_check_square_stochastic( P[[a]] )
		} else {
			 error_msg <- mdp_check_square_stochastic( P[,,a] )
		}
		if (length(error_msg) == 0) {
			a <- a + 1
		} else {
			a <- a1 + 1
		}
	}
}

# Check of R
if (!is_error_detected)  {
	if (is.list(R)) {
		s3 <- dim(R[[1]])[1]
		s4 <- dim(R[[1]])[2]
		a2 <- length(R)
	} else if (inherits(R, 'sparseMatrix-class')) {
		s3 <- dim(R)[2]
		s4 <- s3
		a2 <- dim(R)[1]
	} else if (length(dim(R)) == 3) {
		s3 <- dim(R)[1]
		s4 <- dim(R)[2]
		a2 <- dim(R)[3]
	} else {
		s3 <- dim(R)[1]
		a2 <- dim(R)[2]
		s4 <- s3
	}
	if (s3 < 1 | a2 < 1 | s3 != s4) {
		error_msg <- 'MDP Toolbox ERROR: The reward matrix R must be an array (S,S,A) or (SxA) with S : number of states greater than 0 and A : number of actions greater than 0';
		is_error_detected <- T
	}
}

if (!is_error_detected) {
	if (s1 != s3 | a1 != a2) {
		error_msg <- 'MDP Toolbox ERROR: Incompatibility between P and R dimensions'
		is_error_detected <- T
	}
}

return(error_msg)

}

