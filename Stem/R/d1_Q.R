`d1_Q` <-
function(n,X,d1_X,sigma2omega,B){          
	d1_Q <- n*sum(diag(solve(X) %*% d1_X))-(1/sigma2omega)*sum(diag(solve(X) %*% d1_X %*% solve(X) %*% B))
	return(d1_Q)
}

