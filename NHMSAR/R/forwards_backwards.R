forwards_backwards <-
function(
    prior, transmat, obslik, filter_only = FALSE
){
	T = dim(obslik)[2]
	M = length(prior)
	
	scale = matrix(1,1,T)
	loglik = 0
	alpha = matrix(0,M,T)
	gamma = matrix(0,M,T)
	xi = array(0,c(M,M,T-1))
	t = 1 # correspond to t=order+1
	alpha[,1] = prior * obslik[,t]
	scale[t] = sum(alpha[,t]) 
	alpha[,t]  = normalise(alpha[,t])
	
	transmat2 = t(transmat)
	
	# FORWARD ..............................................................
	for (t in 2:T) {
		tmp = (transmat2 %*% alpha[,t-1]) * obslik[,t]
		scale[t] = sum(tmp) 
		alpha[,t] = normalise(tmp) 
		if (filter_only){ 
			xi[,,t-1] = normalise((alpha[,t-1] %*% t(obslik[,t])) * transmat) }
		}
		loglik = sum(log(scale))
	if (filter_only) { 
		gamma = alpha
	} else {
		beta = matrix(0,M,T)
		gamma = matrix(0,M,T)
		beta[,T] = matrix(1,M,1)
		gamma[,T] = normalise(alpha[,T] * beta[,T])
		t=T
		
		# BACKWARD .........................................................
		for (t in seq(T-1,1,-1)) {
			b = beta[,t+1] * obslik[,t+1];
			beta[,t] = normalise((transmat %*% b))
			gamma[,t] = normalise(alpha[,t] * beta[,t])
			xi[,,t] = normalise((transmat * (alpha[,t] %*% t(b))))
		}
	}
    list(gamma = gamma, xi = xi, loglik = loglik,M=M,alpha=alpha,beta=beta)
}
