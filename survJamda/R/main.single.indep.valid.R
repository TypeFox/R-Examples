main.single.indep.valid <-
function(geno.files, surv.data,normalization = "zscore", method = "none",gn.nb = 100, perf.eval = "auc")
{
	if (!is.element(normalization, c("zscore","combat")))
		stop("\rnormalization = \"zscore\" or normalization = \"combat\"", call.=FALSE)
	
	if(normalization == "combat")
		batchID = det.batchID(geno.files)
	
	for (i in 1:length(geno.files)){
		for (j in 1:length(geno.files))
			if (i != j){	
				common.gene = intersect(colnames(get(geno.files[i])), colnames(get(geno.files[j])))
				ds1 = excl.missing.single.indep(geno.files,i,surv.data,common.gene)
				ds2 = excl.missing.single.indep(geno.files,j,surv.data,common.gene)
		
				if (normalization == "combat")
					mat = prepcombat.single.indep(ds1$mat,ds2$mat,i,j,batchID)
				else
					mat = prepzscore(ds1$mat,ds2$mat)

				
				i.adj = mat[1:nrow(ds1$mat),]
				j.adj = mat[(nrow(ds1$mat)+1):nrow(mat),]
				cat ("Train data set: ", geno.files[j], " Test data set: ", geno.files[i], "\n")
				calPerformance.single.indep(list(mat=j.adj,phyno=ds2$phyno),list(mat=i.adj, phyno=ds1$phyno), method=method,gn.nb=gn.nb, perf.eval)
			}
	}
}

