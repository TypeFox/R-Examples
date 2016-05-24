`d12_Q` <-
function(n,X,d1_X_theta,d1_X_logb,sigma2omega,B){   
	d12_Q <- -n*sum(diag(solve(X) %*% d1_X_theta %*% solve(X) %*% d1_X_logb)) + (2/sigma2omega)*sum(diag(solve(X) %*% d1_X_theta
             %*% solve(X) %*% d1_X_logb %*%  solve(X) %*% B))
	return(d12_Q)
}

