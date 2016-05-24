`d2_Q` <-
function(n,X,d1_X,d2_X,sigma2omega,B){    
	d2_Q <- n*sum(diag(solve(X) %*% d2_X))-n*sum(diag(solve(X) %*% d1_X %*% solve(X) %*% d1_X))-(1/sigma2omega)*sum(diag(solve(X) %*% d2_X %*% solve(X) %*% B))+(2/sigma2omega)*sum(diag(solve(X) %*% d1_X %*% solve(X) %*% d1_X %*% solve(X) %*% B))
	return(d2_Q)
}

