group.graph <-
function(V,cv=NULL,p1,p2=NULL)
{
#only develop for undirected graph
	adjm=diag(0,V)
	allp=t(combn(1:V,2))
	ce=rbinom(n=choose(V,2),size=1,prob=p1)
	adjm[allp]=ce
	if(length(cv)!=0)
	{
		gallp=t(combn(cv,2))
		gce=rbinom(n=choose(length(cv),2),size=1,prob=p2)
		adjm[gallp]=0
		adjm[gallp]=gce
	}
	ind <- lower.tri(adjm) 
	adjm[ind] <- t(adjm)[ind] 
	g=graph.adjacency(adjm, mode=c("undirected"))
	return(g)
}
