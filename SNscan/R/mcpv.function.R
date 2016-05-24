mcpv.function <-
function(obs.stat,ms.stat,direction)
{
	D=get(direction)
	pv=apply(t(t(obs.stat)),1,function(x) (sum(D(ms.stat,x))+1)/(length(ms.stat)+1))
	return(pv=pv)
}
