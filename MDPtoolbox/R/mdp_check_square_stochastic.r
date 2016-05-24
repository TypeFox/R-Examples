mdp_check_square_stochastic <- function(X) {

error_msg <- ''

s1 <- dim(X)[1]
s2 <- dim(X)[2]

if (s1 != s2) {
	error_msg <- 'MDP Toolbox ERROR: Matrix must be square'
} else if ( max(abs(rowSums(X) - rep(1,s2))) > 10^(-12) ) {
	error_msg <- 'MDP Toolbox ERROR: Row sums of the matrix must be 1'
} else if (length(which(X < 0)) > 0) {
	error_msg <- 'MDP Toolbox ERROR: Probabilities must be non-negative'
}

return(error_msg)

}

# Note the difference between Matlab and R: rowSums(Z) in R is sum(Z,2) in Matlab
