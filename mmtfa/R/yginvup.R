yginvup <-
function(p,G,yg){
	yginv <- array(0, dim=c(p, p, G))
	for(g in 1:G){
		diag(yginv[,,g]) <- diag(yg[,,g]^(-1))
	}
	yginv
}
