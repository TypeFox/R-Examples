norm.stat <-
function(obs,pop=1,zloc)
{
	N=length(obs)
	m0=sum(obs)/N;v0=var(obs)
	lnL0= -N*log(sqrt(2*pi))-N*log(sqrt(v0))-sum((obs-m0)^2/(2*v0))
	mz=sum(obs[zloc])/length(obs[zloc]);lz=sum(obs[-zloc])/length(obs[-zloc])
	vz=(sum((obs[zloc]-mz)^2)+sum((obs[-zloc]-lz)^2))/N
	lnL1= -N*log(sqrt(2*pi))-N*log(sqrt(vz))-N/2
	t.stat=N*log(sqrt(v0))+sum((obs-m0)^2/(2*v0))-N*log(sqrt(vz))-N/2
	#equivalent to minimize vz
	return(c(t.stat,m0,mz))
}
