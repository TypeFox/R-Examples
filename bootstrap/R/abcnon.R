"abcnon" <-
function(x, tt,  epsilon = 0.001, alpha = c(.025,.05,.1,.16,.84,.9,.95,.975))
{
call <- match.call()
#abc confidence intervals for nonparametric problems

#tt(P ,x) is statistic in resampling form, where P[i] is weight on x[i]

        if(is.matrix(x)) {n <- nrow(x)} else {n <- length(x)}

	ep <- epsilon/n; I<- diag(n); P0<- rep(1/n,n)
	t0 <- tt(P0,x)	
#calculate t. and t..  .................................................
	t. <- t.. <- numeric(n)
	for(i in 1:n) {  di <- I[i,  ] - P0
			 tp <- tt(P0 + ep * di,x)
			 tm <- tt(P0 - ep * di,x)
			 t.[i] <- (tp - tm)/(2 * ep)
			 t..[i] <- (tp - 2 * t0 + tm)/ep^2}
#calculate sighat,a,z0,and cq ..........................................
	sighat <- sqrt(sum(t.^2))/n
	a <- (sum(t.^3))/(6 * n^3 * sighat^3)	
	delta <- t./(n^2 * sighat)
	cq <- (tt(P0+ep*delta,x) -2*t0 + tt(P0-ep*delta,x))/(2*sighat*ep^2)
	bhat <- sum(t..)/(2 * n^2)
	curv <- bhat/sighat - cq
	z0 <- qnorm(2 * pnorm(a) * pnorm( - curv))	
#calculate interval endpoints............................................
	Z <- z0 + qnorm(alpha)
	za <- Z/(1 - a * Z)^2
	stan <- t0 + sighat * qnorm(alpha)
	abc <- seq(alpha)
        pp <- matrix(0,nrow=n,ncol=length(alpha))
	for(i in seq(alpha)) {abc[i] <- tt(P0 + za[i] * delta,x)
                              pp[,i] <- P0 + za[i] * delta
       }
	limits <- cbind(alpha, abc, stan)
        dimnames(limits)[[2]] <- c("alpha", "abc", "stan")
#output in list form.....................................................
        return(list(limits=limits, 
                    stats=list(t0=t0,sighat=sighat,bhat=bhat), 
                    constants=list(a=a,z0=z0,cq=cq), 
                    tt.inf=t., 
                    pp=pp, 
                    call=call))
}
