null.par = function(K, y, X) 
{
    y = y - mean(y) 
	A = function(lambda)  lambda*K + diag(1, length(y))
    obj = function(lambda)  -length(y)*log(s2.e(lambda)) - log(det(A(lambda)))         

    if (missing(X)) s2.e = function(lambda)  mean(y*solve(A(lambda), y))
    else {
    	mean.X = apply(X, 2, mean)
    	X = sweep(X, 2, mean.X) 
    	b = function(lambda) { XA = t(solve(A(lambda), X)); solve(XA %*% X, XA %*% y) }
    	s2.e = function(lambda) { res = y - X %*% b(lambda); mean(res*solve(A(lambda), res)) }
    }
    	
    tmp = optimize(obj, interval=c(0,1), maximum=TRUE)
    lambda = tmp$maximum

    if (missing(X)) list(y1=y, Sigma0.1 = solve(s2.e(lambda)*A(lambda)))
    else list(y1 = as.vector(y - X %*% b(lambda)), Sigma0.1 = solve(s2.e(lambda)*A(lambda)))
}


