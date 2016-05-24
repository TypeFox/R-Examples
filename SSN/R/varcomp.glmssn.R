varcomp <-
function(object)
{
	r2fe <- GR2(object)
	cp <- covparms(object)
	vcs <- c(r2fe, (1 - r2fe)*cp[cp[,"Parameter"] == "parsill","Estimate"]/
		sum(cp[cp[,"Parameter"] == "parsill","Estimate"]))
	vcnames <- c("Covariates (R-sq)", as.character(cp[cp[,"Parameter"] == 
		"parsill","Covariance.Model"]))
	return(data.frame(VarComp = vcnames, Proportion = vcs))
}


