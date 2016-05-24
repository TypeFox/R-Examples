`ReadSymMatrixFromTriangle` <-
function(file, n.vars){

	cov.vec <- c(scan(file))

	temp <- matrix(0,n.vars,n.vars)
		
	temp[upper.tri(temp,diag=TRUE)] <- cov.vec  # Yes, upper!

	diagonal <- diag(temp)

	resulting.matrix <- temp + t(temp)
	
	diag(resulting.matrix) <- diagonal

	resulting.matrix

} # closes the function

