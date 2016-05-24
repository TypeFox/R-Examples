
splitZscore2.merge.indep <-
function (common.gene, geno.files,surv.data,lst, i,j,method,gn.nb,perf.eval, normalization)
{
	train.ind = NULL

	cat ("Train data sets: ")
	for (m in i)
		cat (geno.files[m], " ")

	train.ind = 1:nrow(lst$mat)

	test.phyno = comb.surv.censor(geno.files,c(j),surv.data)

	test = get(geno.files[j])[!is.na(test.phyno$surv),common.gene]
	test = scale(t(scale(t(test))))

	lst$mat = rbind (lst$mat, test)

	test.phyno$censor = test.phyno$censor[!is.na(test.phyno$surv)]
	test.phyno$surv = test.phyno$surv[!is.na(test.phyno$surv)]

	lst$phyno$surv = c(lst$phyno$surv, test.phyno$surv)
	lst$phyno$censor = c(lst$phyno$censor, test.phyno$censor)

	cat("Test data set: ", geno.files[j], "\n")
	
	test.ind = det.set.ind(geno.files,0,j)


	calPerformance.merge.indep(lst, train.ind, test.ind, method,gn.nb,perf.eval, normalization)
}

