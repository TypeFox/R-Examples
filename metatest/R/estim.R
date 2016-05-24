estim <-
function(y, x, v, threshold=0.00001, maxiter=100) {
	
	# ML estimation derived form Viechtbauer code
	# RETURNS: sigmasqn, the between study variance, the number of iterations,
	# and the convergence flag
	change <- 1.0
	iter <- 0
	sigmasq <- 0
	ifault <- 0    # becomes one if max iter is exceeded, becomes -1 if var estimate is set to zero
	
	sigmasq	<- var(y) - mean(v)
	
	intc <- dim(x)[2]
	
	while (change > threshold) {
		iter=     iter + 1
		ww=       (1/(v+sigmasq))
		invsumww= 1/sum(ww)
		w=        diag(ww)				
		wsq=      w^2
		# 	normal model
		if(intc>0) {
			b=        (solve(t(x) %*% w %*% x)) %*% t(x) %*% w %*% y
			hlp=      t(((y-x%*%b)^2 )-v) %*% wsq
		} else { # model without variables
			b=matrix(0,1,1)
			hlp=      t(((y)^2 )-v) %*% wsq
		}
		
		# b =  (solve(t(x) %*% w %*% x)) %*% t(x) %*% w %*% y
		
		sigmasqn= (sum(hlp)/sum(diag(wsq)) ) # ML
		change=   abs(sigmasqn - sigmasq)
		sigmasq=  sigmasqn
		if (sigmasqn  < 0) {sigmasqn=0;ifault=-1;break}
		if (iter > maxiter) {ifault=1;break}
	}
	
	out=c(sigmasqn,iter,ifault)
	
	return(out)
}

