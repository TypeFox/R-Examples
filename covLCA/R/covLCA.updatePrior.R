covLCA.updatePrior <-
function(b,xx,R) 
{
    bb <- matrix(b,ncol=(R-1)) #A: rows=covariates, col=LC : beta_pj
    exb <- exp(xx %*% bb) #A: rows=indiv, col=LC
    p <- cbind(exb,1)/(rowSums(exb)+1) #A: LC probabilities, rows=indiv, cols=LC (where last column corresponds to last LC, the reference) 
	
	return(p)
}
