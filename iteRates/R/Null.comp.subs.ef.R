Null.comp.subs.ef<-function(treelist,mod.id=c(1,0,0,0),trace){
	loop<-length(treelist)
	FPR <- NA
	res<-list()
	for (i in 1:loop){
		if(trace){if(i%%200==0){cat("\n",i,"of",loop)}
		if(i%%10==0){cat("*")}}
		res[[i]]<-comp.subs(treelist[[i]],mod.id=mod.id,verbose=FALSE)[,c(14,19)]
		FPR[i] <- sum(res[[i]][,1]<.05,na.rm=T)/sum(!is.na(res[[i]][,1]))	
		}
	return(list(FPR=FPR,res=res))
	}
