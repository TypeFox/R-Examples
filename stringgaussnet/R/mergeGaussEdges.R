mergeGaussEdges <-
function (Edges1,Edges2,Names)
{
	Common1<-rep(FALSE,nrow(Edges1))
	Common2<-rep(FALSE,nrow(Edges2))
	for (i in 1:nrow(Edges1))
	{
		node1<-Edges1$node1[i]
		node2<-Edges1$node2[i]
		CommonPos<-which((Edges2$node1==node1 & Edges2$node2==node2) | (Edges2$node2==node1 & Edges2$node1==node2))
		if (length(CommonPos)>0){Common1[i]<-TRUE;Common2[CommonPos]<-TRUE}
	}
	
	Edges1M<-Edges1[,c("node1","node2","Rho","P.Value")]
	Edges1M$Network<-rep(NA,nrow(Edges1M))
	Edges1M$Network[Common1]<-"Common"
	Edges1M$Network[!Common1]<-Names[1]
	
	Edges2M<-Edges2[!Common2,c("node1","node2","Rho","P.Value")]
	Edges2M$Network<-rep(Names[2],nrow(Edges2M))
	
	MergedEdges<-rbind(Edges1M,Edges2M)
	return(MergedEdges)
}
