padj.qpcr <-
function(data,controls=NULL,method="BH") {
	for (n in controls) {
		nn=grep(n,row.names(data))
		if (sum(grep("pval.z",names(data)))>0) data$pval.z[nn]=NA
		data$pval.mcmc[nn]=NA
	}
	if (sum(grep("pval.z",names(data)))>0) data$padj.z=p.adjust(data$pval.z,method=method)
	data$padj.mcmc=p.adjust(data$pval.mcmc,method=method)
	return(data)
}
