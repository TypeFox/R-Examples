histPheno <-
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
		

		layout.mat <- get.layout.mat(dim(data.obj$pheno)[2])
		layout(layout.mat)
		for(p in 1:length(pheno.labels)){
			hist(data.obj$pheno[,p], xlab = pheno.labels[p], main = pheno.labels[p])
			}
	}
