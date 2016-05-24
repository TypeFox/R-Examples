Likelihood_thetas <-
function(thetas,covmat) {
		cholcov <- chol(covmat)
		invcholcov <- ginv(cholcov)
		logsqrtdet <- sum(log(diag(cholcov)))
		-(1/2)*(colSums((crossprod(invcholcov,thetas))^2))-logsqrtdet
	}
