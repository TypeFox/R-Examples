snps.to.genes <-
function(snp.info,gene.info,distance){
	j=1
	gene.snps=c()
	for(j in 1:nrow(gene.info)){
		snps.inds=snp.info$Chr==gene.info$Chr[j] & snp.info$Position>=(gene.info[j,"Start"]-distance) & snp.info$Position<=(gene.info[j,"End"]+distance)
		gene.snps=c(gene.snps,list(snp.info$Name[which(snps.inds=="TRUE")]))
	}
	rm(j)
	names(gene.snps)=gene.info$Name
	gene.snps
	}
