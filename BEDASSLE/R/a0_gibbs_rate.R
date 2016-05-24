a0_gibbs_rate <-
function(thetas,covmat,a0){
		cholcov <- chol(covmat)
		invcholcov <- ginv(cholcov)
		tmp <- (1/2)*(colSums((crossprod(invcholcov,thetas))^2))
		return(sum(tmp)/a0)
	}
