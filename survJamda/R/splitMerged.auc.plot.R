splitMerged.auc.plot <-
function (geno.files,lst, i,j,col, method, time.dep)
{
	normalization = ifelse(col == "black", "ComBat", "Zscore1")
	cat ("Normalization = ", normalization, "\n")

	cat ("Train data sets: ")
	train.ind = det.set.ind(geno.files,1,i)

	cat("Test data set: ", geno.files[j], "\n")

	test.ind = det.set.ind(geno.files,0,j)

	calPerformance.auc.plot(lst, train.ind, test.ind, geno.files[j],col, method, normalization, time.dep)
}

