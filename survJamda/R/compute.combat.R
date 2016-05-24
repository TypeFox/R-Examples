compute.combat <-
function (fileGeno, fileSample)
{
	mat.adj  = ComBat(fileGeno, fileSample, write = FALSE, filter = FALSE, prior.plots= FALSE, skip = 1, par.prior = TRUE)
	gn.names = rownames(mat.adj)
	mat.adj = as.matrix(mat.adj[,2:ncol(mat.adj)])
	rownames(mat.adj) = gn.names

	return (mat.adj)
}

