compareSubTrees <-
function(trees,focal,k,mod.id=c(1,0,0,0),min.val=0.01){#compares trees for a given k
	#trees should be a list of all subtrees
	#focal should be a vector of subtrees to compare
	if(is.list(trees)==FALSE)stop("\nSubtrees must be object class list")
	if(length(trees)<2)stop("\nMultiple subtrees required")
	tempmin<-NA
		subTax<-1:length(focal)#letters[1:length(focal)]#here's the id letters for the subtrees
		treeComp<-c(trees[focal])#list of subtrees to compare
			no.trees<-length(treeComp)
				for(Tr in 1:no.trees){#
					
					thisTre<-treeComp[[Tr]]
					tempmin[Tr]<-max(branching.times(thisTre))*min.val
					treeComp[[Tr]]$b<-thisTre$edge[,2]>length(thisTre$tip.label)
					treeComp[[Tr]]$y<-thisTre$edge.length
					}
			min.branch<-min(tempmin)		
		PosClus<-posClusters(subTax,k)
		compcluS<-list()
			for (PC in 1:length(PosClus)){
				compcluS[[PC]]<-calcClus(PosClus[[PC]],treeComp,min.branch=min.branch,mod.id=mod.id)
				}
		return(compcluS)
		}

