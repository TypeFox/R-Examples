Psi.lnorm <- function(x, p, eta, sq.var, theta.fix, theta.var)
{
	res <- 0
	for(i in 1:length(eta))
		for(j in 1:length(eta))
			if(p[i,j] != 0)
				res <- res + p[i,j] * psi.lnorm(x, eta[[i]], eta[[j]], sq.var[[i]], sq.var[[j]], theta.fix[[i]], theta.var[[i,j]])
	res
}