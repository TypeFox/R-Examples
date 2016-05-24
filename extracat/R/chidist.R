

chidist <- function(x,dim = 1){
	
	# sum profile
	s <- apply(x,-dim,sum)
	
	sx <- apply(x,dim,function(z) z/sum(z))
	sx <- sx/sqrt(s)
	
	D <- dist(t(sx))
	return(D)
	
}

chidist2 <- function(M){
	cs <- colSums(M)
	
	M2 <- M/rowSums(M)
	
	S <- M2 %*% diag(1/cs) %*% t(M2)
	
	s <- diag(S)
	
	D <- s%*% t(rep(1,length(s))) + rep(1,length(s)) %*% t(s) - 2*S
	return(D)
}



# equality

# coords <- ca1$rowcoord %*% diag((ca1$sv))
# dist(coords)/sqrt(N)