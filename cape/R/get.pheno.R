get.pheno <-
function(data.obj, scan.what = c("eigentraits", "raw.traits")){
	
		#If the user does not specify a scan.what, 
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentraits, otherwise, use raw phenotypes
	type.choice <- c(grep("eig", scan.what), grep("ET", scan.what), grep("et", scan.what)) #look for any version of eigen or eigentrait, the user might use.
	if(length(type.choice) > 0){ #if we find any, use the eigentrait matrix
		pheno <- data.obj$ET
		if(is.null(pheno)){stop("There are no eigentraits. Please set scan.what to raw.traits, or run get.eigentraits().")}
		}else{
			pheno <- data.obj$pheno #otherwise, use the raw phenotype matrix
			}
	
	return(pheno)
	
	
}
