## compute congruence coefficients (or other FUN) over the last dimension for a 3-dim array
## returns an array with one more row and column and the same number of layers

congruence.array <-
	function(X, FUN=congruence.coef, stat.name="phi", round=FALSE, scale=1, ref="last"){
	
	if (length(dim(X)) != 3) stop("Not a 3-way array")
	d1 <- dim(X)[1]
	d2 <- dim(X)[2]
	d3 <- dim(X)[3]
	names <- dimnames(X)
	
	base <- switch(ref, first=1, last=d3, ref)
	if (!is.numeric(base) && base <= d3) stop("ref level out of range")
	 	

# TODO: no need to separate code for d3 == 2, >2
	if (d3 == 2) {
		
		Y <- array(NA, c(d1+1, d2+1, d3))
		Y[1:d1,1:d2,1:d3] <- X
		
		for (i in 1:d1) { Y[i,d2+1,] <- FUN(X[i,,1], X[i,,2]) }
		for (j in 1:d2) { Y[d1+1,j,] <- FUN(X[,j,1], X[,j,2]) }
		Y[d1+1,d2+1,] <- FUN(as.vector(X[,,1]),as.vector(X[,,2]))
		
		}
	
	if (d3 >= 3) {
		
		Y <- array(NA, c(d1+1, d2+1, d3))
		Y[1:d1,1:d2,1:d3] <- X
		
		for (i in 1:d1){
			for (k in 1:d3){
				Y[i,d2+1,k] <- FUN(X[i,,k], X[i,,base])
				}
			}
		
		for (j in 1:d2){
			for (k in 1:d3){
				Y[d1+1,j,k] <- FUN(X[,j,k], X[,j,base])
				}
			}

		for (k in 1:d3){
			Y[d1+1,d2+1,k] <- FUN(as.vector(X[,,k]),as.vector(X[,,base]))
			}

		}
		
	# assign dimnames
	dimnames(Y) <- list(
		c(names[[1]], stat.name),
		c(names[[2]], stat.name),
		names[[3]]
		)
	Y <- scale * Y

	if (is.logical(round) && round) round(Y)
	else if (is.numeric(round)) round(Y,digits=round)
	else Y
	
	}