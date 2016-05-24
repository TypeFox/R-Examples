get.eigentraits <-
function(data.obj, scale.pheno = TRUE, normalize.pheno = TRUE){

	#first make sure there are no individuals
	#with missing phenotypes. This also makes 
	#sure the phenotypes are numeric
	
	pheno <- data.obj$pheno
	ind.names <- rownames(pheno)
	
	#make sure all the phenotypes are numeric
	pheno <- apply(pheno, 2, as.numeric)
	rownames(pheno) <- ind.names
	
	#see if there are any individuals with missing 
	#phenotypes data. Remove individuals with 
	#missing data from the pheno matrix.

	na.rows <- which(is.na(pheno), arr.ind = TRUE)
	rows.to.remove <- sort(unique(na.rows[,1]))

	#stop and warn if all individuals are about to be removed
	if(length(rows.to.remove) == dim(pheno)[1]){
		stop("All individuals have missing phenotype data. Try narrowing your phenotype selection with select.pheno()")
		}
	
	if(length(rows.to.remove) > 0){	
		message("Removing ", length(rows.to.remove), " individuals due to missing phenotypes.\n")
		data.obj <- remove.ind(data.obj, ind.which = rows.to.remove)
		}

	pheno <- data.obj$pheno

	#if there are fewer phenotypes than individuals...
	if(dim(pheno)[1] < dim(pheno)[2]){
			message("\nThere must be more individuals than phenotypes.\nPlease select fewer phenotypes to analyze or\nperform a dimension reduction on the\nphenotypes before reading in the data.")
			return(data.obj)
			}

	#This function mean centers and standardizes a vector
	center.std <- function(v){
		mean.v <- mean(v)
		centered <- v - mean.v
		sd.v <- sd(v)
		final.v <- centered/sd.v
		return(final.v)
		}

	if(normalize.pheno){
		pheno <- apply(pheno, 2, rz.transform)
		}

	if(scale.pheno){
		pheno <- apply(pheno, 2, center.std) #mean center and standardize the phenotypes
		}


	data.obj$pheno <- pheno #replace the raw phenotypes with scaled, normalized phenotypes (if we have done those things)

	svd.pheno <- svd(pheno)
	

	#add the eigentraits and singular values to the data object
	ET <- svd.pheno$u
	rownames(ET) <- rownames(pheno); colnames(ET) <- paste("ET", 1:dim(ET)[2], sep = "")
	data.obj$ET <- ET
	data.obj$singular.values <- svd.pheno$d
	data.obj$right.singular.vectors <- svd.pheno$v
	

	
	return(data.obj)
	
}
