network.mc.scan <-
function(n,g,radius=3,attribute,model,pattern,fix.edge=FALSE,max.prop=0.5,xmin=NULL,zetatable=NULL)
{
	if (pattern!="structure")
	{
		#permute attribute
		stat = get(model, mode = "function")
		if(model=="norm.stat"|model=="multinom.stat"| model=="powerlaw.stat"| model=="conpowerlaw.stat")
		{
			attsim.data=NULL
			for(k in 1:n)
			{
				s=sample(attribute$obs,length(attribute$obs))
				attsim.data=rbind(attsim.data,s)
			}
		} else
		{
			p=sum(attribute$obs)/sum(attribute$pop)
			attsim.data=matrix(rbinom(n=n*length(attribute$pop), 
				size=attribute$pop, prob=p),byrow=TRUE,nrow=n)
		}
	}
	if (pattern!="attribute")
	{
		#permute structure
		Sg=graph.rmedge(n=n,g=g,fix.edge)
	}
	maxk=NULL
	for(s in 1:n)
	{
		s.attribute=attribute
		if (pattern!="structure") s.attribute$obs=attsim.data[s,] else s.attribute$obs=attribute$obs
		if (pattern!="attribute") sg=Sg[[s]] else sg=g
		mc=network.scan(g=sg,radius=radius,attribute=s.attribute,model=model,pattern=pattern,max.prop,xmin,zetatable)
		if(length(nrow(mc))==0) maxk=rbind(maxk,mc) else maxk=rbind(maxk,mc[1,])
	}
	return(maxk)
}
