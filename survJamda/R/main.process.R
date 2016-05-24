main.process <-
function(common.gene, geno.files, surv.data,method = "none", time.dep)
{
	curr_set = 1:length(geno.files)
	norm.set = c("combat","zscore1", "zscore2")
	color.set = c("black","blue","green")
	train.lst = NULL

	batchID = det.batchID(geno.files)
	
	for (y in curr_set){
		x = setdiff(curr_set, y)
		for (normalization in norm.set){
			prep = get(paste("prep",normalization, sep = ""))
			if(normalization == "combat")
				lst = prep(common.gene,geno.files,surv.data,batchID, x,y)
			else
				lst = prep(common.gene,geno.files,surv.data, x,y)	
			switch (normalization,
				"combat"=(col="black"),
				"zscore1"=(col = "blue"),
				"zscore2"=(col = "green")
			)

			if (normalization == "zscore1" || normalization == "combat")
				splitMerged.auc.plot (geno.files, lst, x, y,col, method, time.dep)
			else
				splitZscore2.auc.plot (common.gene,geno.files,surv.data, lst, x, y,col, method, time.dep)
		}
		if (time.dep) abline(0.5,0,col = "red") 
		if (time.dep) legend(50,0.2,legend = c("ComBat","Zscore1","Zscore2"), text.col = color.set, bty = "n", pt.cex = .5) 
	}

	if (time.dep) title(main = "Time-dependent AUC based on predicted time\nIndependent validation",outer = TRUE, cex.main=1.5) 
}

