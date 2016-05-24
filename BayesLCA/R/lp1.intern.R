lp1.intern <-
function(parvec, dat, counts.n, rm.ind, prior.list, val.ind, G0, M0){

	alpha<- prior.list$alpha
	beta<- prior.list$beta
	delta<- prior.list$delta

	theta<- rep(0, G0*M0)
	len.parvec.theta<- G0*M0 - length(rm.ind)

	if(length(rm.ind)>0) theta[-rm.ind]<- parvec[1:len.parvec.theta] else theta<- parvec[1:len.parvec.theta]
	theta[rm.ind]<- val.ind

	theta<- matrix(theta, G0, M0)
	tau<- parvec[(len.parvec.theta+1):(length(parvec))]

	lpsum<- 0
	
	for(g in 1:(G0-1)) lpsum<- lpsum +  tau[g]*apply(X=dat, MARGIN=1, FUN=fx.g.intern, theta1=theta[g,])
	
	lpsum<- lpsum +  (1-sum(tau))*apply(X=dat, MARGIN=1, FUN=fx.g.intern, theta[G0,])
	sum(counts.n*log(lpsum) + log(ddirichlet(c(tau, 1-sum(tau)), delta)) + sum(dbeta(theta, alpha, beta, log=TRUE)) )
}
