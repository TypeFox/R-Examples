get.nb<-function(arc,pal,index=NULL)
{
	if(is.null(index))
		index<-1:length(pal[[1]][[1]])

	lindex<-length(index)
	
	nb<-vector(mode="list", length=lindex)

	arcid<-arc[[1]]$ArcId
	lpoly<-arc[[1]]$LeftPoly
	rpoly<-arc[[1]]$RightPoly

	for(p in 1:lindex)
	{
		arcs<-sort(unique(abs(pal[[2]][[ index[p] ]][[1]])))

		if(arcs[1]==0)
			arcs<-arcs[-1]

		arcindex<-as.vector(tapply( arcs , 1:length(arcs),
			function(X,arcid){which(arcid==X)}
			,arcid) )

		thisnb<-c(lpoly[arcindex], rpoly[arcindex])
		nb[[p]]<-sort(unique(thisnb))
	}

	nb
}
