pred.time.indep.valid <-
function(geno.files, surv.data, method = "none", time.dep = 0)
{
	common.gene = colnames(get(geno.files[1]))
	for (i in 2:length(geno.files))
		common.gene = intersect(common.gene, colnames(get(geno.files[i])))

	par (mfrow = c(1,length(geno.files)), cex = 1)
	par(oma=c(.5,.5,length(geno.files),.5))

	main.process (common.gene, geno.files, surv.data,  method = "none", time.dep)
}

