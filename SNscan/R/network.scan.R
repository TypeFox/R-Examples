network.scan <-
function(g,radius=3,attribute,model,pattern,max.prop=0.5,xmin=NULL,zetatable=NULL)
{
	#determine the shortest pathes among nodes
	sp=shortest.paths(g)
	#construct scanning regions: splist[[center]][[path]] 
	splist=NULL
	for(d in 1:radius)
	{
		spi=apply(sp,1,function(x) list(which(x<=d)))
		if(d==1) splist=spi else splist=mapply(splist, spi, FUN=c, SIMPLIFY=FALSE)
	}
	#only attribute clustering
	if (pattern=="attribute")
	{
		stat = get(model, mode = "function")
		obs.Mcluster = NULL
		for (k in 1:length(splist)) 
		{
			if(model=="powerlaw.stat"|model=="conpowerlaw.stat")inoutm=cbind(rep(k,length(splist[[k]])),1:length(splist[[k]]),
				t(sapply(splist[[k]],function(x) c(stat(obs=attribute$obs,pop=attribute$pop,zloc=x,xmin=xmin,zetatable=zetatable),length(x)))))
			else inoutm=cbind(rep(k,length(splist[[k]])),1:length(splist[[k]]),
				t(sapply(splist[[k]],function(x) c(stat(obs=attribute$obs,pop=attribute$pop,zloc=x),length(x)))))
			obs.Mcluster=rbind(obs.Mcluster,inoutm)
			if(model!="multinom.stat") colnames(obs.Mcluster)=c("C","D","test.L","m0","mz","z.length")
			else colnames(obs.Mcluster)=c("C","D","test.L",paste("P0",unique(attribute$obs),sep=""),
					paste("P1",unique(attribute$obs),sep=""),"z.length")
			obs.Mcluster[is.nan(obs.Mcluster)]=0
			if(nrow(obs.Mcluster)==1)obs.rankinfo=obs.Mcluster else {
				obs.rankinfo=obs.Mcluster[order(obs.Mcluster[,3],decreasing=TRUE),]}
			drow=which(obs.rankinfo[,ncol(obs.rankinfo)]>floor(length(V(g))*max.prop)|obs.rankinfo[,ncol(obs.rankinfo)]==1)			
			if (length(drow)==0) obs.rankinfo=obs.rankinfo else obs.rankinfo=obs.rankinfo[-drow,]
		}
	#only structure clustering
	} else if (pattern=="structure"){
		obs.Mcluster = NULL
		for (k in 1:length(splist)) 
		{
			inoutm=cbind(rep(k,length(splist[[k]])),1:length(splist[[k]]),
				t(sapply(splist[[k]],function(x) c(structure.stat(g,subnodes=x),length(x)))))
			obs.Mcluster=rbind(obs.Mcluster,inoutm)
			colnames(obs.Mcluster)=c("C","D","test.L","S0","Sz","z.length")
			obs.Mcluster[is.nan(obs.Mcluster)]=0
			obs.rankinfo=obs.Mcluster[order(obs.Mcluster[,3],decreasing=TRUE),]
			drow=which(obs.rankinfo[,ncol(obs.rankinfo)]>floor(length(V(g))*max.prop)|
				obs.rankinfo[,ncol(obs.rankinfo)]==1|obs.rankinfo[,4]>obs.rankinfo[,5])			
			if (length(drow)==0) obs.rankinfo=obs.rankinfo else obs.rankinfo=obs.rankinfo[-drow,]
		}
	#for both clustering	
	} else if (pattern=="both"){
		stat = get(model, mode = "function")
		obs.Mcluster = NULL
		for (k in 1:length(splist)) 
		{
			inoutm=cbind(rep(k,length(splist[[k]])),1:length(splist[[k]]),
				t(sapply(splist[[k]],
				function(x) c(structure.stat(g,subnodes=x)[1]+stat(obs=attribute$obs,pop=attribute$pop,zloc=x)[1],
				structure.stat(g,subnodes=x),stat(obs=attribute$obs,pop=attribute$pop,zloc=x),
				length(x)))))
			obs.Mcluster=rbind(obs.Mcluster,inoutm)
			if(model!="multinom.stat") colnames(obs.Mcluster)=c("C","D","test.L","test.S","S0","Sz",
				"test.A","m0","mz","z.length")
			else colnames(obs.Mcluster)=c("C","D","test.L","test.S","S0","Sz",
					"test.A",paste("P0",unique(attribute$obs),sep=""),
					paste("P1",unique(attribute$obs),sep=""),"z.length")
			obs.Mcluster[is.nan(obs.Mcluster)]=0
			obs.rankinfo=obs.Mcluster[order(obs.Mcluster[,3],decreasing=TRUE),]
			drow=which(obs.rankinfo[,ncol(obs.rankinfo)]>floor(length(V(g))*max.prop)|obs.rankinfo[,ncol(obs.rankinfo)]==1)			
			if (length(drow)==0) obs.rankinfo=obs.rankinfo else obs.rankinfo=obs.rankinfo[-drow,]
		}
	}
	return(obs.rankinfo)
}
