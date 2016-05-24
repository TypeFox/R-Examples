padjustOTU=function(otu.summary,method="BH") {
	sigs=c()
	otus=names(otu.summary$otuWise)
	sn=otu.summary
	pvs=c()
	for(i in otus) {
		tub=sn$otuWise[[i]]
		ntr=length(dimnames(tub)$pvalue)
		for(th in ntr:2) {
			for(tv in 1:(th-1)) {
				pvs=append(pvs,tub[th,tv])
			}
		}
	}
	pvs=p.adjust(pvs,method=method)
	counts=0
	for(i in otus) {
		tub=sn$otuWise[[i]]
		ntr=length(dimnames(tub)$pvalue)
		for(th in ntr:2) {
			for(tv in 1:(th-1)) {
				counts=counts+1
				sn$otuWise[[i]][th,tv]=pvs[counts]
			}
		}
	}
	return(sn)
}