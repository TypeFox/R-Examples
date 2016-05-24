CompareSubTrees.k <-
function(trees,focal,maxk,mod.id=c(1,0,0,0),min.val=0.01){	res<-list()
	Ks<-1:maxk
		for (i in 1:maxk){
			cat("\n",i,"of k =",maxk)
			res[[i]]<-compareSubTrees(trees,focal=focal,k=Ks[i],mod.id=mod.id,min.val=min.val)
			}
			cat("\n")
		return(res)
		}

