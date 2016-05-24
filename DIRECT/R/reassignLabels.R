reassignLabels <-
function (x, pars)
{
	labels = unique (x)
	nlabels = length (labels)
	x.new = x
	pars.new = matrix (0, nrow=nlabels, ncol=ncol(pars))

	for (i in 1:nlabels)
	{
		x.new[which (x==labels[i])] = i
		pars.new[i,] = pars[labels[i],]
	}
	
	return (list (c=x.new, pars=pars.new))
}

