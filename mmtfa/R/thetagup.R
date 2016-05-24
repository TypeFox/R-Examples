thetagup <-
function(q,G,betag,lg,sg,thetag){
	for(g in 1:G){
		if(q==1){
			thetag[,,g] <- diag(q) - betag[,,g] %*% lg[,,g] + 
					betag[,,g] %*% sg[,,g] %*% betag[,,g]
		}
		if(q>1){
			thetag[,,g] <- diag(q) - betag[,,g] %*% lg[,,g] + 
					betag[,,g] %*% sg[,,g] %*% t(betag[,,g])
		}
	}
	thetag
}
