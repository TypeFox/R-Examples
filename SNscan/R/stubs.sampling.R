stubs.sampling <-
function(s,g)#g: permuting graph
{
	k=degree(g)
	v=c(1:length(k))
	vs=rep(v,k)
	S=NULL
	for(i in 1:s)
	{
		S=rbind(S,sample(vs,length(vs)))
	}
	return(S)
}
