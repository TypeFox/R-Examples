relabel <-
function (probs.mcmc, nIter, nItem, nClust, RELABEL.THRESHOLD, PRINT=0, PACKAGE="DIRECT") {
	probs.mcmc = as.matrix (probs.mcmc)
	probs = as.double (as.vector (t(probs.mcmc)))
	result = .C ("relabel_R", perm_mtx = as.integer (rep (0, length=nClust*nIter)), probs.mcmc, nIter, as.integer(nItem), as.integer(nClust), as.double (RELABEL.THRESHOLD), as.integer(PRINT), PACKAGE=PACKAGE)
	perms.mtx = matrix (result$perm_mtx, byrow=TRUE, ncol=nClust)
	return (perms.mtx)
}

