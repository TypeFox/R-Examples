sigmaup <-
function(p,G,lg,yg,sigma){
	for(g in 1:G){
		sigma[,,g] <- lg[,,g] %*% t(lg[,,g]) + yg[,,g]
	}
	sigma
}
