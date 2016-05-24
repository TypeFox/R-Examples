evaluate.weights.2.lnorm <- function(x, w, p, eta, sq.var, theta.fix, theta.var, weights.evaluation.epsilon, weights.evaluation.max.iter, support.epsilon)
{
	w <- rep(1, length(x)) / length(x)
	
	KLp.curr <- 1
	KLp.prev <- 0
	
	iter <- 0
	while( (KLp.curr / KLp.prev > 1 + weights.evaluation.epsilon) && (iter < weights.evaluation.max.iter) )
	{
		iter <- iter + 1
		delta_ <- 0
		for(i in 1:length(eta))
			for(j in 1:length(eta))
				if(p[i,j] != 0)
				{
					ev.eta.i <- eta[[i]](x, theta.fix[[i]]) 
                    sq.sigma.i <- log(1 + sq.var[[i]](x, theta.fix[[i]]) / (ev.eta.i * ev.eta.i))
                    mu.i <- log(ev.eta.i) - 0.5 * sq.sigma.i

					theta.var[[i,j]] <- optim(
						par = theta.var[[i,j]],
						function(theta) KLD.new(x, w, ev.eta.i, sq.sigma.i, mu.i, eta[[j]], sq.var[[j]], theta)
						)$par
					
					delta_ <- delta_ + p[i,j] * psi.lnorm(x, eta[[i]], eta[[j]], sq.var[[i]], sq.var[[j]], theta.fix[[i]], theta.var[[i,j]])
				}
			
		k.max <- which.max(delta_)
		k.min <- which.min(delta_)

		KLp.prev <- sum(w * delta_)

		op.res <- optimize(
			f = function(alpha) 
			{
				delta_ <- 0
				w.t <- w_(alpha, w, k.max, k.min)
				g(x, w.t, p, eta, sq.var, theta.fix, theta.var)
			},
			interval = c(0, w[k.min]),
			maximum = TRUE)
		w <- w_(op.res$maximum, w, k.max, k.min)
		KLp.curr <- op.res$objective
		x <- x[w > support.epsilon]
		w <- w[w > support.epsilon]
	}
	list(x = x, w = w, theta.var = theta.var)
}