beta.pd.utils <- function(com, tree) {
	combin	<-	combn(nrow(com),2)
	labcomb <-	apply(combin,2,function(x) paste(rownames(com)[x],collapse="-"))
	pd.obs	<-	pd2(com,tree)[,"PD"]
	com.tot <- t(apply(combin,2,function(x) colSums(com[x,])>0))
	pd.obs.tot <- pd2(com.tot,tree)[,"PD"]
	sum.pd.obs <- apply(combin,2,function(x) sum(pd.obs[x]))
	min.pd.obs <- apply(pd.obs.tot-t(combn(pd.obs,2)),1,min)
	dif.pd.obs <- apply(combin,2,function(x) diff(pd.obs[x]))
	return(list(pd.obs=pd.obs,pd.obs.tot=pd.obs.tot,sum.pd.obs=sum.pd.obs,min.pd.obs=min.pd.obs,dif.pd.obs=dif.pd.obs,labcomb=labcomb,combin=combin))
}

