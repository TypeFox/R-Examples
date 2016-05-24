sigmainvup <-
function(p,G,yginv,lg,q,sigmainv){
	for(g in 1:G){
		sigmainv[,,g] <- yginv[,,g] - yginv[,,g] %*% lg[,,g] %*% solve(diag(q) + 
				t(lg[,,g]) %*% yginv[,,g] %*% lg[,,g]) %*% t(lg[,,g]) %*% yginv[,,g]

	}
	sigmainv
}
