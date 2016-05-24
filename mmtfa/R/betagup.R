betagup <-
function(q,p,G,lg,sigmainv,betag){
	for(g in 1:G){
		betag[,,g] <- t(lg[,,g]) %*% sigmainv[,,g]
	}
	betag
}
