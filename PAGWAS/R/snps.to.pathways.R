snps.to.pathways <-
function(pathways,gene.snps){
	  snps.paths.list=as.list(rep(NA,length(pathways)))
		for(p in 1:length(pathways)){
			snps=c()
			snps=unique(unlist(gene.snps[names(gene.snps) %in% pathways[[p]]]))
			snps.paths.list[[p]]=snps
		}
		names(snps.paths.list)=names(pathways)
		snps.paths.list
}
