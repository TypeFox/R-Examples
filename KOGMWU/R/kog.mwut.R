kog.mwut <-
function(gos,Alternative="t") {
	terms=unique(gos$term)
	gos$seq=as.character(gos$seq)
	nrg=gos[!duplicated(gos$seq),2]
	names(nrg)=gos[!duplicated(gos$seq),1]
	rnk=rank(nrg)
	names(rnk)=names(nrg)
	pvals=c();drs=c();nams=c();levs=c();nseqs=c()
	for (t in terms){
		got=gos[gos$term==t,]
		got=got[!duplicated(got$seq),]
		ngot=gos[gos$term!=t,]
		ngot=ngot[!duplicated(ngot$seq),]
		ngot=ngot[!(ngot$seq %in% got$seq),]
		sgo.yes=got$seq
		n1=length(sgo.yes)
		sgo.no=ngot$seq
		n2=length(sgo.no)
		if (n2 < n1) {
			print(paste("skipping",t,"nseqs =",n1))
			next
		}
		wi=wilcox.test(nrg[sgo.yes],nrg[sgo.no],alternative=Alternative)	 
		r1=sum(rnk[sgo.yes])/n1
		r0=sum(rnk[sgo.no])/n2
		dr=r1-r0
		drs=append(drs,round(dr,0))
		levs=append(levs,got$lev[1])
		pvals=append(pvals,wi$p.value)
		nseqs=append(nseqs,n1)	
	}
	res=data.frame(cbind(nseqs,"delta.rank"=drs,"pval"=pvals))
	res=cbind("term"=as.character(terms),res)
	res$pval=as.numeric(as.character(res$pval))
	res$delta.rank=as.numeric(as.character(res$delta.rank))
	res$nseqs=as.numeric(as.character(res$nseqs))
	res=res[order(res$pval),]
	res$padj=p.adjust(res$pval,method="BH")
	return(res)
}
