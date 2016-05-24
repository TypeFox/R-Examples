generate.survival.data <-
function(gene.nb, tot.genes, sample.nb, beta.init, correlation, shape, scale){

	my_beta = rep(0,tot.genes)
	my_beta[1:(gene.nb/2)] = beta.init
	my_beta[((gene.nb/2)+1):gene.nb]=-beta.init
	
	u = runif(sample.nb)

 	X = corgen(len = sample.nb,r = correlation)
	ds1 = X$x

	for (i in 1:(tot.genes-1))
		 ds1 = cbind(ds1,corgen(x = X$x, r = correlation)$y)

	lambda = exp(as.vector(my_beta%*%t(ds1)))
	div = (-log (1-u))/(scale*lambda)
	T = div^(1/shape)

	v = runif(n=length(T),min=min(T), max = quantile(T, probs = seq(0, 1, 0.9))[2])
	 censor =  as.numeric (T <= v)
  	 T = pmin(T,v)
	
	cat("\nSingle data set\n")
	iter.crossval(ds1, T, censor,ngroup = 10, gn.nb = gene.nb,zscore =0,gn.nb.display = 0)
	return(list(ds1=ds1,T=T,censor=censor))
}

