CreateFormula <-
function(data, i, covariable=NULL) {
	if(length(covariable)==0) { fmla <- as.formula(data[,1] ~ data[,i]); return(fmla)
	} else {
		covname <- colnames(covariable)
		fmla <- as.formula(paste("data[,1] ~ data[,i] + ", paste(covname, collapse= "+")))
	}	
	return(fmla)
}
