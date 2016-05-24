discrandz <-
function(n,G){
	zmat <- matrix(0, n, G)
	for(i in 1:n){
		zmat[i, sample(1:G, 1)] <- 1
	}
	zmat
}
