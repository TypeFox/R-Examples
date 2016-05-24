qqPheno <-
function(data.obj, pheno.which = NULL, pheno.labels = NULL){
	
	all.pheno <- data.obj$pheno


	if(is.null(pheno.which)){
		pheno.names <- colnames(all.pheno)
		}else{
			if(is.numeric(pheno.which)[1]){
				pheno.names <- colnames(all.pheno)[pheno.which]
				}else{
				pheno.names <- pheno.which	
				}
			}
	
	if(is.null(pheno.labels)){
		pheno.labels <- pheno.names
		}
		
	pheno.num.pairs <- pair.matrix(1:length(pheno.labels))

		layout.mat <- get.layout.mat(dim(pheno.num.pairs)[1])
		layout(layout.mat)
		for(p in 1:nrow(pheno.num.pairs)){
			qqplot(data.obj$pheno[,pheno.num.pairs[p,1]], data.obj$pheno[,pheno.num.pairs[p,2]], xlab = pheno.labels[pheno.num.pairs[p,1]], ylab = pheno.labels[pheno.num.pairs[p,2]])
			}
	}
