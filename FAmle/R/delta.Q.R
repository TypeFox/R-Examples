delta.Q <-
function(p,model,ln=FALSE)
{
	k <- model$k
	if(k > 1)
	{
		cov.hat <- model$cov.hat
		prime <- sapply(as.list(1:k),function(h) Diff.3(h,model,p,ln))
		VAR <- diag(cov.hat)
		COV <- cov.hat[lower.tri(cov.hat)]
		Q.var.hat <- sum(2*COV*sapply(as.list(as.data.frame(combn(1:k,2))),function(g) prod(prime[g]))) +
			sum(VAR*prime^2)
		if(!ln) return(c(mu.hat=distr(p,model=model,type='q'),sd.hat=sqrt(Q.var.hat)))
		else return(c(mu.hat=log(distr(p,model=model,type='q')),sd.hat=sqrt(Q.var.hat)))
	}
	else
	{
		h <- 1e-4
		if(!ln)
		{
			Q.var.hat <- model$cov.hat[1,1]*((distr(x=p,dist=model$dist,param=model$par.hat+h,type='q')-
				distr(x=p,dist=model$dist,param=model$par.hat-h,type='q'))/2/h)^2
			Q.hat <- distr(x=p,model=model,type='q')
		}
		else
		{
			Q.var.hat <- model$cov.hat[1,1]*((log(distr(x=p,dist=model$dist,param=model$par.hat+h,type='q'))-
				log(distr(x=p,dist=model$dist,param=model$par.hat-h,type='q')))/2/h)^2
			Q.hat <- log(distr(x=p,model=model,type='q'))
		}
		return(c(mu.hat=Q.hat,sd.hat=sqrt(Q.var.hat)))
	}
}
