norm.pheno <-
function(data.obj, mean.center = TRUE){
	
	
	pheno <- data.obj$pheno
	raw.pheno <- data.obj$pheno #retain the raw phenotype values

	#This function mean centers and standardizes a vector
	center.std <- function(v){
		mean.v <- mean(v, na.rm = TRUE)
		centered <- v - mean.v
		sd.v <- sd(v, na.rm = TRUE)
		final.v <- centered/sd.v
		return(final.v)
		}

	pheno <- apply(pheno, 2, rz.transform)


	if(mean.center){
		pheno <- apply(pheno, 2, center.std) #mean center and standardize the phenotypes
		}



	data.obj$pheno <- pheno #replace the raw phenotypes with scaled, normalized phenotypes (if we have done those things)
	data.obj$raw.pheno <- raw.pheno

	return(data.obj)
	
	
	
}
