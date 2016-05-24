binom.stat <-
function(obs,pop,zloc)
{
	p0=sum(obs)/sum(pop);p1=sum(obs[zloc])/sum(pop[zloc])
	p10=(sum(obs)-sum(obs[zloc]))/(sum(pop)-sum(pop[zloc]))
	t.stat=(sum(obs[zloc])*log(p1/p0)+(sum(pop[zloc])-sum(obs[zloc]))*log((1-p1)/(1-p0))+
		sum(obs[-zloc])*log(p10/p0)+(sum(pop[-zloc])-sum(obs[-zloc]))*log((1-p10)/(1-p0)))*(p1>p10)+
		1*(p1<=p10)
	#equivalent to minimize vz
	return(c(t.stat,p10,p1))
}
