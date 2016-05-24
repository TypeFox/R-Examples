emtest.thm3 <-
function(tb,N=10000,tol=1e-8)
#This function computes a_h values in Theorem 3 of Li and Chen (2009) 
#
#tb:  covariance matrix \tilde B_{22}.
#N:   number of repetitions.
#tol: a value below which is judged as 0.
#
#Need to first install the R package quadprog.  
{	
	m0=dim(tb)[1]	
	eig.tb=eigen(tb)
	tb2=eig.tb$vectors%*%diag(sqrt(eig.tb$values))%*%t(eig.tb$vectors)

	output=c()
	for(i in 1:N)
	{
		out=solve.QP(Dmat=tb,dvec=tb2%*%rnorm(m0,0,1),Amat=diag(rep(1,m0)),bvec=rep(0,m0))
		output=c(output,sum(out$solution>tol))
	}
	ah=output
	table(ah)/N
}
