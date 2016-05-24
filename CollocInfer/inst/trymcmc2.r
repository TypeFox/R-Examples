source('Devel/Colloc.MCMC.R')

walk.var = diag(c(0.001,0.001,0.01))
cscale = 0.000001

prior = function(pars){ return( sum(dnorm(pars/c(1,1,10),log=TRUE))) }
nstep = 100


hmm = Colloc.MCMC(times,data,pars+c(-0.1,0.1,-1),coefs,lik,proc,prior,walk.var,cscale,nstep,in.meth='SplineEst',control.in=control.in)