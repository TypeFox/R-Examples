main.merge.indep.valid <-
function(geno.files,surv.data,gn.nb=100,method = "none", normalization= "zscore1",perf.eval = "auc")
{
	if (length(geno.files) < 3) 
		stop ("\rThere should be minimum 3 data sets", call. = FALSE)

	if (!is.element(normalization, c("combat","zscore1","zscore2")))
		stop ("\rSet normalization = \"combat\" or normalization = \"zscore1\" or normalization = \"zscore2\"", call. = FALSE)

	if (normalization == "combat")
		batchID = det.batchID(geno.files)

	common.gene = colnames(get(geno.files[1]))
	for (i in 2:length(geno.files))
		common.gene = intersect(common.gene, colnames(get(geno.files[i])))

	curr_set = 1:length(geno.files)

	for (y in curr_set){
		x = setdiff(curr_set, y)
		prep = get(paste("prep",normalization, sep = ""))

		if (normalization == "combat")
			lst = prep(common.gene,geno.files,surv.data,batchID,x,y)		
		else
			lst = prep(common.gene,geno.files,surv.data,x,y)

		if (normalization == "zscore1" || normalization == "combat")
			splitMerged.indep (geno.files,lst, x, y, method,gn.nb,perf.eval, normalization)
					
		else
			splitZscore2.merge.indep (common.gene,geno.files,surv.data,lst, x, y, method,gn.nb,perf.eval, normalization)	
	}
}

