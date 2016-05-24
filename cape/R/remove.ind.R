remove.ind <-
function(data.obj, ind.which){
	
	#take out the missing entries from the phenotype object
	data.obj$pheno <- data.obj$pheno[-ind.which,,drop=FALSE]

	if(!is.null(data.obj$p.covar.table)){
		data.obj$p.covar.table <- data.obj$p.covar.table[-ind.which,,drop=FALSE]
		}
	if(!is.null(data.obj$g.covar.table)){
		data.obj$g.covar.table <- data.obj$g.covar.table[-ind.which,,drop=FALSE]
		}

	#check the dimensions of the phenotype object and remove
	#individuals accordingly
	if(!is.null(data.obj$raw.pheno)){
		data.obj$raw.pheno <- data.obj$raw.pheno[-ind.which,]
		}
	#remove individuals from the genotype matrix if necessary
	if(!is.null(data.obj$geno)){
		data.obj$geno <- data.obj$geno[-ind.which,]
		}

	
	return(data.obj)
	
}
