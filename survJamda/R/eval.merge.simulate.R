eval.merge.simulate <-
function(d1,d2,tot.genes, gene.nb, zscore)
{
	mat = rbind(d1$ds1,d2$ds1)
	
	cat("\nMerged data set\n")
	iter.crossval(mat, c(d1$T,d2$T), c(d1$censor,d2$censor),ngroup = 10,zscore =1,gn.nb =gene.nb,gn.nb.display = 0)
}

