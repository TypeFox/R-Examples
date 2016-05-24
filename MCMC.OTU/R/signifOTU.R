signifOTU=function(otu.summary,p.cutoff=0.05) {
# otu.summary=qss;p.cutoff=0.05
	sigs=c()
	otus=names(otu.summary$otuWise)
	sn=otu.summary
	for(i in otus) {
		tub=sn$otuWise[[i]]
		ntr=length(dimnames(tub)$pvalue)
		for(th in ntr:2) {
			for(tv in 1:(th-1)) {
				if (tub[th,tv]<=p.cutoff) { 
					if (!(i %in% sigs))	{sigs=append(sigs,i)}
				}
			}
		}
	}
	return(sigs)
}