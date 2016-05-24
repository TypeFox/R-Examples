delete.pheno <-
function(data.obj, phenotypes){
	
	pheno <- data.obj$pheno
	
	pheno.col <- get.col.num(pheno, phenotypes)
	
	if(length(pheno.col) > 0){
		new.pheno <- pheno[,-pheno.col]
		data.obj$pheno <- new.pheno
		}
	
	return(data.obj)
	
}
