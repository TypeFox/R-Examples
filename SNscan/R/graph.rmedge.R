graph.rmedge <-
function(n,g,fix.edge=TRUE)
{
	k=degree(g)
	nV=length(V(g))
	allp=t(combn(nV,2))
	pe=k[allp[,1]]*k[allp[,2]]/(sum(k)) #upper triangular probability
	adjm=diag(0,nV)
	Sg=NULL
	for(i in 1:n)
	{
		if (fix.edge==TRUE)ce=rmulti.one(size=sum(k)/2,p=pe) else ce=rbinom(n=nrow(allp),size=1,prob=pe)
		adjm[allp]=ce
		ind <- lower.tri(adjm) 
		adjm[ind] <- t(adjm)[ind] 
		sg=graph.adjacency(adjm, mode=c("undirected"))
		Sg=c(Sg,list(sg))
	}
	return(Sg)
}
