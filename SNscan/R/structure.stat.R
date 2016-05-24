structure.stat <-
function(g,subnodes)#g: igraph object 
{
	subg=induced.subgraph(g, subnodes)
	CZ=sum(degree(subg))/2;CG=sum(degree(g))/2
	MuZ=(sum(degree(g)[subnodes]))^2/(4*CG);MuG=CG
	t.stat=CZ*log(CZ/MuZ)+(CG-CZ)*log((CG-CZ)/(MuG-MuZ))
	p1=(CZ/MuZ);p10=(CG-CZ)/(MuG-MuZ)
	return(c(t.stat,p10,p1))
}
