det.set.meta <-
function (i,j, geno.files,surv.data, method)
{
	cat ("Train data sets: ")
	for (m in i)
		cat (geno.files[m], " ")
		
	cat("Test data set: ", geno.files[j], "\n")
	common.gene = colnames(get(geno.files[i[1]]))

	for (k in 2:length(i))
		common.gene = intersect(common.gene, colnames(get(geno.files[i[k]])))

	zstat = pool.zscores(common.gene, i, geno.files,surv.data)
	calPerformance.meta(common.gene,zstat, i, j, geno.files,surv.data,method)
}

