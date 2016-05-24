create.pathway.df <-function(genotypes,snps.paths){
	 mat=matrix(0,nrow=ncol(genotypes),ncol=length(snps.paths))
	 rownames(mat)=colnames(genotypes)
	 colnames(mat)=names(snps.paths)
	 for(p in 1:length(snps.paths)){
		inds=which(rownames(mat) %in% snps.paths[[p]])
		mat[inds,p]=1
	 }
	 data.frame(mat)
	}
