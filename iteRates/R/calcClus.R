calcClus <-
function(PosCl,treeComp,min.branch,mod.id){

	res<-list()
	reS<-array(dim=c(length(PosCl$group[,1]),6))
	k<-length(PosCl$VecBin)
	textGroup<-GroupsAsText(PosCl)
		TESTS<-list()
		yS<-list()
		bS<-list()
		starts<-c(cumsum(c(1,PosCl$VecBin)))
		starts<-starts[-length(starts)]
		ends<-cumsum(PosCl$VecBin)
		for (i in 1:PosCl$dimen[2]){
			yS[[i]]<-treeComp[[i]]$y
			bS[[i]]<-treeComp[[i]]$b
			}
			
		for (i in 1:PosCl$dimen[1]){
			k <- length(PosCl$VecBin)
			groups <- paste(textGroup[i,],collapse=" vs. ")
				for (j in 1:length(PosCl$VecBin)){
				y<-unlist(yS[PosCl$group[i,starts[j]:ends[j]]])
				b<-unlist(bS[PosCl$group[i,starts[j]:ends[j]]])
				TESTS[[j]]<-model.sel.distbn(b,y,min.branch=min.branch,mod.id=mod.id)
				}
			res[[i]]<-list(k=k,groups=groups,tests=TESTS)
		}
		return(res)
}

