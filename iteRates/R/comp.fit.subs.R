comp.fit.subs <-
function(trees,focal,k,mod.id=c(1,0,0,0),min.val=0.01){
		cst<-CompareSubTrees.k(trees=trees,focal=focal,maxk=k,mod.id=mod.id,min.val=min.val)
		tt<-tab.table(cst)
		res<-tt.Summary(tabtabres=tt,trees=trees,focal=focal)
		return(res)
	}

