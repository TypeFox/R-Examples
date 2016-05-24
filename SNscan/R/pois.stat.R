pois.stat <-
function(obs,pop,zloc)
{
	p0=sum(obs)/sum(pop);p=sum(obs[zloc])/sum(pop[zloc]);Q=sum(obs[-zloc])/sum(pop[-zloc])
	t.stat=(sum(obs[zloc])*log(p/p0)+sum(obs[-zloc])*log(Q/p0))*(p>Q)+1*(p<=Q)
	#equivalent to minimize vz
	return(c(t.stat,Q,p))
}
