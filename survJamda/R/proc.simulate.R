proc.simulate <-
function(tot.genes = 100, correlation = 0, gene.nb = 50, sample.nb = 400, beta.init = 0.5, shape = 1, scale = 1)
{
	d1 = generate.survival.data (gene.nb, tot.genes,sample.nb,beta.init, correlation, shape, scale)
	old.beta.init = beta.init
	old.shape = shape
	old.scale = scale

	cat ("Enter a different value for beta.init for the second data set. Enter only a numeric value. Example: 2\n")
	beta.init = as.numeric(readline())
	if (!is.numeric(beta.init))
		beta.init = old.beta.init
	
	cat ("Enter a different value for shape for the second data set. Enter only a numeric value. Example: 2\n")
	shape = as.numeric(readline())
	if (!is.numeric(shape))
		shape = old.shape

	cat ("Enter a different value for scale for the second data set. Enter only a numeric value. Example: 1.5\n")
	scale = as.numeric(readline())
	if (!is.numeric(scale))
		scale = old.scale

	d2 = generate.survival.data(gene.nb, tot.genes,sample.nb,beta.init, correlation, shape, scale)

	eval.merge.simulate(d1,d2,tot.genes, gene.nb, zscore = 1)
}

