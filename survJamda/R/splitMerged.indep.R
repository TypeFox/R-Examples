splitMerged.indep <-
function (geno.files,lst, i,j, method,gn.nb,perf.eval, normalization)
{
	cat ("Train data sets: ", geno.files[i], " ")
	train.ind = lst$train.ind

	cat("Test data set: ", geno.files[j], "\n")

	test.ind = (length(lst$train.ind)+1):nrow(lst$mat)

	calPerformance.merge.indep(lst, train.ind, test.ind, method,gn.nb,perf.eval, normalization)
}

