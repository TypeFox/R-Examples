contrandz <-
function(n,G){
	zmat <- matrix(0, n, G)
	for(i in 1:n){
		for(g in 1:G){
			zmat[i, g] <- runif(1, min=0, max=1)
		}
	}
	zmat <- zmat/rowSums(zmat)
	zmat
}
