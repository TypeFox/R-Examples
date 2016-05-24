LASSO.exact = function(lambda, tau, X, y)
{
 	n = nrow(X)
 	p = ncol(X)
    X = apply(X, 2, scale)
    X = X/sqrt(n-1)
    y = scale(y)/sqrt(n-1)
    lambda = lambda - 1E-10
 
    beta = rep(0, p)
    breaks = rep(0, p)
    next.tau = tau
    change = rep("|", p)
 	
 	if(all(abs(tau)<0.5)){            # no selected variables
 		S = t(X)%*%y
 		breaks = S 
 		next.tau[order(S, decreasing=TRUE)[1]] = 1
 		change[order(S, decreasing=TRUE)[1]] = "+"
 	} else {
 		ins = (abs(tau) > 0.5)
 		X1 = as.matrix(X[, ins])
 		X2 = as.matrix(X[, !ins])
 		sign1 = tau[abs(tau) > 0.5]
 
 		tmp = solve(t(X1)%*%X1, cbind(t(X1)%*%y, sign1))
 		beta1.int = tmp[,1]
 		beta1.b = -tmp[,2] 
 		ols.res = y - X1 %*% beta1.int
 		S2.int = t(X2)%*%ols.res
 		S2.b = t(X2)%*% (X1%*%tmp[,2])
     	s1 = S2.int/(sign(S2.b)-S2.b)
     	s2 = -S2.int/(sign(S2.b)+S2.b)

     	breaks[ins] = -beta1.int/beta1.b            
 		breaks[!ins] = abs(sign(S2.b)*ifelse(s1>s2 & s1<lambda, s1, -s2))  
        next.lambda = max(c(0, breaks[breaks<lambda]))
        C1 = !(breaks<next.lambda) & breaks<lambda
    	beta[ins] = beta1.int + beta1.b*next.lambda
        S = t(X)%*%(y-X1%*%beta[ins])                   # score at next lambda
        next.tau[C1 & ins] = 0
        next.tau[C1 & !ins] = sign(S[C1])
        change[C1] = "+"
        change[C1 & ins] = "-"
 	}

 	data.frame(beta = beta, score = S, breaks = breaks, tau = next.tau, change = change)
}
  
