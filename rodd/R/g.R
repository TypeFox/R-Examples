g <- function(x, w, p, eta, sq.var, theta.fix, theta.var)
{
	res <- 0
	for(i in 1:length(eta))
		for(j in 1:length(eta))
			if(p[i,j] != 0)
			{
				ev.eta.i <- eta[[i]](x, theta.fix[[i]]) 
                        sq.sigma.i <- log(1 + sq.var[[i]](x, theta.fix[[i]]) / (ev.eta.i * ev.eta.i))
                        mu.i <- log(ev.eta.i) - 0.5 * sq.sigma.i

				opt.res <- optim(
					par = theta.var[[i,j]],
					function(theta) KLD.new(x, w, ev.eta.i, sq.sigma.i, mu.i, eta[[j]], sq.var[[j]], theta)
					)
				theta.var[[i,j]] <- opt.res$par
				res <- res + p[i,j] * opt.res$value
			}
	res
}