kog.ft <-
function(gos) {
	terms=unique(gos$term)
	gos$seq=as.character(gos$seq)
	pft=c();nam=c();lev=c();nseqs=c()
	for (t in terms) {
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		go.sig=sum(got$value)
		go.ns=length(got[,1])-go.sig
		ngo.sig=sum(ngot$value)
		ngo.ns=length(ngot[,1])-ngo.sig
		sig=c(go.sig,ngo.sig) # number of significant genes belonging and not belonging to the tested GO category
		ns=c(go.ns,ngo.ns) # number of not-significant genes belonging and not belonging to the tested GO category
		mm=matrix(c(sig,ns),nrow=2,dimnames=list(ns=c("go","notgo"),sig=c("go","notgo")))
		ff=fisher.test(mm,alternative="greater")
		pft=append(pft,ff$p.value)
		nam=append(nam,as.character(got$name[1]))
		lev=append(lev,got$lev[1])
		nseqs=append(nseqs,length(got[,1]))
	}
	res=data.frame(cbind("term"=as.character(terms),nseqs,"pval"=pft))
	res$pval=as.numeric(as.character(res$pval))
	res$nseqs=as.numeric(as.character(res$nseqs))
	res=res[order(res$pval),]
	res$padj=p.adjust(res$pval,method="BH")
	return(res)
}
